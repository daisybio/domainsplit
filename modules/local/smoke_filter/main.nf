process SMOKE_FILTER {
    tag "smoke_filter"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path cobinet_db
    val  n_ddis

    output:
    path "cobinet.smoke.sqlite3", emit: cobinet_db
    path "versions.yml",          emit: versions

    script:
    """
    #!/usr/bin/env python3
    import shutil
    import sqlite3
    import sys

    N = ${n_ddis}
    HALF = N // 2

    shutil.copy("${cobinet_db}", "cobinet.smoke.sqlite3")
    con = sqlite3.connect("cobinet.smoke.sqlite3")
    con.execute("PRAGMA foreign_keys=ON")

    pos_before = con.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction WHERE negative=0"
    ).fetchone()[0]
    neg_before = con.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction WHERE negative=1"
    ).fetchone()[0]
    dom_before = con.execute("SELECT COUNT(*) FROM domain").fetchone()[0]
    print(
        f"smoke filter: before -> pos={pos_before} neg={neg_before} domain={dom_before}; "
        f"keep {HALF} of each",
        flush=True,
    )

    con.execute("CREATE TEMP TABLE keep_ddi(id INTEGER PRIMARY KEY)")
    con.execute(
        "INSERT INTO keep_ddi(id) "
        "SELECT id FROM domain_domain_interaction WHERE negative=0 "
        "ORDER BY RANDOM() LIMIT ?",
        (HALF,),
    )
    con.execute(
        "INSERT INTO keep_ddi(id) "
        "SELECT id FROM domain_domain_interaction WHERE negative=1 "
        "ORDER BY RANDOM() LIMIT ?",
        (HALF,),
    )
    con.execute(
        "DELETE FROM domain_domain_interaction "
        "WHERE id NOT IN (SELECT id FROM keep_ddi)"
    )
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
        f"smoke filter: after  -> pos={pos_after} neg={neg_after} domain={dom_after}",
        flush=True,
    )
    con.close()

    con = sqlite3.connect("cobinet.smoke.sqlite3")
    con.execute("VACUUM")
    con.close()

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\\n")
    """
}
