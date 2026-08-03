process FILTER_DB {
    tag "filter_db"
    label 'process_high'
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
    import time

    shutil.copy("${domainsplit_db}", "domainsplit.filter.sqlite3")
    con = sqlite3.connect("domainsplit.filter.sqlite3")
    con.execute("PRAGMA foreign_keys=ON")
    con.execute("PRAGMA synchronous=OFF")
    con.execute("PRAGMA journal_mode=MEMORY")
    con.execute("PRAGMA temp_store=MEMORY")

    # Logging
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
        map(tuple, ppi_data[["uniprot_id_a", "uniprot_id_b"]].values.tolist())
    ))

    # uniprot_id,pfam_id -> Strings, saver as keys can change ordering
    mapping_data = pd.read_csv("${meta_mapping}")
    mappings = sorted(set(
        map(tuple, mapping_data[["pfam_id", "uniprot_id"]].values.tolist())
    ))


    # Temporary table to filter mapping -> only keep domains and proteins relevant
    con.execute("CREATE TEMP TABLE keep_mapping(pfam TEXT, prot TEXT)")
    con.executemany(
        "INSERT OR IGNORE INTO keep_mapping VALUES (?, ?)",
        [(d, p) for d, p in mappings],
    )
    con.execute("CREATE INDEX idx_keep_mapping ON keep_mapping(pfam, prot)")

    # Temporary table to filter ppi -> only keep ppis we need
    con.execute("CREATE TEMP TABLE keep_ppi(pua TEXT, pub TEXT)")
    con.executemany(
        "INSERT OR IGNORE INTO keep_ppi VALUES (?, ?)",
        [(pa, pb) for pa, pb in ppi_pairs] + [(pb, pa) for pa, pb in ppi_pairs],
    )
    con.execute("CREATE INDEX idx_keep_ppi ON keep_ppi(pua, pub)")
    con.commit()

    t0 = time.time()

    # Filtering procedure
    con.execute('''
        DELETE FROM protein_protein_interaction
        WHERE NOT EXISTS (
            SELECT 1 FROM keep_ppi k
            JOIN protein pa ON pa.uniprot_id = k.pua
            JOIN protein pb ON pb.uniprot_id = k.pub
            WHERE pa.id = protein_protein_interaction.protein_id_a
            AND pb.id = protein_protein_interaction.protein_id_b
        )
    ''')

    con.execute('''
        DELETE FROM domain_protein_map
        WHERE NOT EXISTS (
            SELECT 1 FROM keep_mapping k
            JOIN domain d ON d.pfam_id = k.pfam
            JOIN protein p ON p.uniprot_id = k.prot
            WHERE d.id = domain_protein_map.domain_id
            AND p.id = domain_protein_map.protein_id
        )
    ''')

    con.execute("DROP TABLE keep_mapping")
    con.execute("DROP TABLE keep_ppi")
    con.commit()

    t1 = time.time()
    print(f"Time for PPI and mapping DELETE + cleanup: {t1 - t0:.1f}s", flush=True)


    # Clean up, check domain 
    t2 = time.time()
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
    t3 = time.time()
    print(f"Time for domain and protein DELETE: {t3 - t2:.1f}s", flush=True)

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

    t4 = time.time()
    con = sqlite3.connect("domainsplit.filter.sqlite3")
    con.execute("VACUUM")
    con.close()
    t5 = time.time()
    print(f"Time for VACUUM: {t5 - t4:.1f}s", flush=True)

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
