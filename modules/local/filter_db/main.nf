process FILTER_DB {
    tag "filter_db"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db
    path meta_ppi
    path meta_mapping

    output:
    path "domainsplit.filter.sqlite3", emit: domainsplit_db
    path "versions.yml",              emit: versions

    script:
    """
    #!/usr/bin/env python3
    import shutil
    import sqlite3
    import sys
    import pandas as pd

    shutil.copy("${domainsplit_db}", "domainsplit.filter.sqlite3")
    con = sqlite3.connect("domainsplit.filter.sqlite3")
    con.execute("PRAGMA foreign_keys=ON")

    // Logging
    ppi_before = con.execute("SELECT COUNT(*) FROM protein_protein_interaction").fetchone()[0]
    map_before = con.execute("SELECT COUNT(*) FROM domain_protein_map").fetchone()[0]
    dom_before = con.execute("SELECT COUNT(*) FROM domain").fetchone()[0]
    prot_before = con.execute("SELECT COUNT(*) FROM protein").fetchone()[0]
    print(
        f"filter db: before -> ppi={ppi_before} mapping={map_before} "
        f"domain={dom_before} protein={prot_before}",
        flush=True,
    )

    # Load ppi metadata containing uniprot_id_a,uniprot_id_b,model_entity_id,local_tar_name,has_af_model,protein_id_a,protein_id_b,path
    ppi_data = pd.read_csv("${meta_ppi}")
    ppi_pairs = sorted(set(
        map(tuple, ppi_data[["protein_id_a", "protein_id_b"]].values.tolist())
    ))

    mapping_data = pd.read_csv("${meta_mapping}")
    mappings = sorted(set(
        map(tuple, mapping_data[["domain_id", "protein_id"]].values.tolist())
    ))

    # Filter ppi table for these specific interactions
    # Need to check ordering in both directions

    # Temporary file 
    con.execute("CREATE TEMP TABLE keep_mapping(dom INTEGER, prot INTEGER)")
    con.executemany(
        "INSERT OR IGNORE INTO keep_mapping VALUES (?, ?)",
        [(int(d), int(p)) for d, p in mappings],
    )


    con.execute("CREATE TEMP TABLE keep_ppi(proteina INTEGER, proteinb INTEGER)")
    con.executemany(
        "INSERT OR IGNORE INTO keep_ppi VALUES (?, ?)",
        [(int(pa), int(pb)) for pa, pb in ppi_pairs],
    )

    con.execute("CREATE TEMP TABLE keep_ppi(proteina INTEGER, proteinb INTEGER)")
    con.executemany(
        "INSERT OR IGNORE INTO keep_ppi VALUES (?, ?)",
        [(int(pb), int(pa)) for pa, pb in ppi_pairs],
    )

    con.commit()

    # Also automatically deletes self-interactions
    con.execute(
        "DELETE FROM protein_protein_interaction "
        "WHERE (protein_id_a, protein_id_b) NOT IN (SELECT proteina, proteinb FROM keep_ppi)"
    )

    con.execute(
        "DELETE FROM domain_protein_map "
        "WHERE (protein_id, domain_id) NOT IN (SELECT prot, dom FROM keep_mapping)"
    )

    con.commit()

    # --- Deduplicate reversed-direction rows ---
    # Ensure no interaction is stored multiple times
    # For domain_domain_interaction: unique key is (domain_id_a, domain_id_b, source),
    # so a "reverse duplicate" is (a,b,source) and (b,a,source) both present.
    dupe_ddi = con.execute('''
        SELECT t1.id AS id_keep, t2.id AS id_drop,
            t1.domain_id_a, t1.domain_id_b, t1.source,
            t1.negative AS neg_keep, t2.negative AS neg_drop
        FROM domain_domain_interaction t1
        JOIN domain_domain_interaction t2
        ON t1.domain_id_a = t2.domain_id_b
        AND t1.domain_id_b = t2.domain_id_a
        AND t1.source = t2.source
        WHERE t1.id < t2.id
    ''').fetchall()

    n_conflict = sum(1 for r in dupe_ddi if r[5] != r[6])
    if n_conflict:
        print(f"WARNING: {n_conflict} reversed DDI pairs disagree on 'negative' label", flush=True)
        for r in dupe_ddi:
            if r[5] != r[6]:
                print(f"  conflict: domains {r[2]}/{r[3]} source={r[4]} "
                    f"neg_keep={r[5]} neg_drop={r[6]}", flush=True)

    drop_ids = [r[1] for r in dupe_ddi]
    con.executemany("DELETE FROM domain_domain_interaction WHERE id = ?", [(i,) for i in drop_ids])
    print(f"dedup ddi: removed {len(drop_ids)} reversed-direction duplicates", flush=True)

    # For protein_protein_interaction: unique key is (protein_id_a, protein_id_b), no source column.
    dupe_ppi = con.execute('''
        SELECT t1.protein_id_a, t1.protein_id_b, t2.protein_id_a, t2.protein_id_b,
            t1.score AS score_keep, t2.score AS score_drop
        FROM protein_protein_interaction t1
        JOIN protein_protein_interaction t2
        ON t1.protein_id_a = t2.protein_id_b
        AND t1.protein_id_b = t2.protein_id_a
        WHERE t1.protein_id_a < t1.protein_id_b
    ''').fetchall()

    n_score_conflict = sum(1 for r in dupe_ppi if r[4] != r[5])
    if n_score_conflict:
        print(f"WARNING: {n_score_conflict} reversed PPI pairs disagree on 'score'", flush=True)

    con.executemany(
        "DELETE FROM protein_protein_interaction WHERE protein_id_a = ? AND protein_id_b = ?",
        [(r[2], r[3]) for r in dupe_ppi],  # drop the (b,a) orientation, keep (a,b) where a<b
    )
    con.commit()
    print(f"dedup ppi: removed {len(dupe_ppi)} reversed-direction duplicates", flush=True)

    # Clean up, check domain 
    con.execute('''
       DELETE FROM domain
        WHERE id NOT IN (
            SELECT domain_id_a FROM domain_domain_interaction
            UNION
            SELECT domain_id_b FROM domain_domain_interaction
            UNION
            SELECT domain_id FROM domain_protein_map
        )
    ''')

    con.execute('''
       DELETE FROM protein
        WHERE id NOT IN (
            SELECT protein_id_a FROM protein_protein_interaction
            UNION
            SELECT protein_id_b FROM protein_protein_interaction
            UNION
            SELECT protein_id FROM domain_protein_map
        )
    ''')
    con.commit()

    ppi_after = con.execute("SELECT COUNT(*) FROM protein_protein_interaction").fetchone()[0]
    map_after = con.execute("SELECT COUNT(*) FROM domain_protein_map").fetchone()[0]
    dom_after = con.execute("SELECT COUNT(*) FROM domain").fetchone()[0]
    prot_after = con.execute("SELECT COUNT(*) FROM protein").fetchone()[0]
    print(
        f"filter db: after  -> ppi={ppi_after} mapping={map_after} "
        f"domain={dom_after} protein={prot_after}",
        flush=True,
    )
    con.close()

    con = sqlite3.connect("domainsplit.filter.sqlite3")
    con.execute("VACUUM")
    con.close()

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\\n")
    """

    stub:
    """
    touch domainsplit.filter.sqlite3
    echo '"${task.process}":' > versions.yml
    echo '    stub: "true"' >> versions.yml
    """
}
