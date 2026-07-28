process FILTER_DDI {
    tag "filter_ddi"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db
    path mapping_data

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

    pos_before = con.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction WHERE negative=0"
    ).fetchone()[0]
    neg_before = con.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction WHERE negative=1"
    ).fetchone()[0]
    dom_before = con.execute("SELECT COUNT(*) FROM domain").fetchone()[0]
    print(
        f"filter ddis: before -> pos={pos_before} neg={neg_before} domain={dom_before}; "
        flush=True,
    )

    # Load mapping, containing protein_id,domain_id
    mapping_data = pd.read_csv("${mapping_data}")
    domain_ids = sorted(set(mapping_df["domain_id"].tolist()))

    # Temporary file 
    con.execute("CREATE TEMP TABLE keep_domain(id INTEGER PRIMARY KEY)")
    con.executemany(
        "INSERT OR IGNORE INTO keep_domain(id) VALUES (?)",
        [(int(d),) for d in domain_ids],
    )

    con.execute(
        "DELETE FROM domain_domain_interaction "
        "WHERE domain_id_a NOT IN (SELECT id FROM keep_domain) "
        "   OR domain_id_b NOT IN (SELECT id FROM keep_domain)"
    )

    # Also drop any self-interactions from domain_domain_interaction
    con.execute(
        "DELETE FROM domain_domain_interaction "
        "WHERE domain_id_a = domain_id_b"
    )
    con.commit()

    # Clean up domains no longer referenced by any surviving DDI
    con.execute('''
        DELETE FROM domain
        WHERE id NOT IN (
            SELECT domain_id_a FROM domain_domain_interaction
            UNION
            SELECT domain_id_b FROM domain_domain_interaction
        )
    ''')
    con.commit()

    pos_after = con.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction WHERE negative=0"
    ).fetchone()[0]
    neg_after = con.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction WHERE negative=1"
    ).fetchone()[0]
    dom_after = con.execute("SELECT COUNT(*) FROM domain").fetchone()[0]
    print(
        f"filter ddis: after  -> pos={pos_after} neg={neg_after} domain={dom_after}",
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
