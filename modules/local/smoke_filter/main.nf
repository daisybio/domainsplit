process SMOKE_FILTER {
    tag "smoke_filter"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path sqlite_3did
    path negatome
    val  n_ddis

    output:
    path "3did.smoke.sqlite3", emit: sqlite
    path "negatome.smoke.tsv", emit: negatome
    path "versions.yml",       emit: versions

    script:
    """
    #!/usr/bin/env python3
    import sqlite3
    import random
    import shutil
    import sys

    random.seed(42)
    N = ${n_ddis}

    # ---- 3did SQLite: keep N random positive DDIs; restrict Domain / domain_length to referenced rows
    shutil.copy("${sqlite_3did}", "3did.smoke.sqlite3")
    con = sqlite3.connect("3did.smoke.sqlite3")
    total = con.execute("SELECT COUNT(*) FROM DDI1").fetchone()[0]
    print(f"3did DDI1 rows before: {total}", flush=True)

    keep = min(N, total)
    rows = con.execute(
        "SELECT domain1, domain2 FROM DDI1 ORDER BY RANDOM() LIMIT ?",
        (keep,)
    ).fetchall()

    con.execute("CREATE TEMP TABLE keep_pairs(d1 TEXT, d2 TEXT, PRIMARY KEY(d1, d2))")
    con.executemany("INSERT OR IGNORE INTO keep_pairs VALUES (?, ?)", rows)

    con.execute('''
        DELETE FROM DDI1
        WHERE NOT EXISTS (
            SELECT 1 FROM keep_pairs k
            WHERE k.d1 = DDI1.domain1 AND k.d2 = DDI1.domain2
        )
    ''')

    domains = {d for pair in rows for d in pair}
    con.execute("CREATE TEMP TABLE keep_domains(name TEXT PRIMARY KEY)")
    con.executemany("INSERT OR IGNORE INTO keep_domains VALUES (?)", [(d,) for d in domains])
    con.execute("DELETE FROM Domain WHERE Name NOT IN (SELECT name FROM keep_domains)")
    con.execute("DELETE FROM domain_length WHERE domain NOT IN (SELECT name FROM keep_domains)")
    con.commit()

    after_ddi = con.execute("SELECT COUNT(*) FROM DDI1").fetchone()[0]
    after_dom = con.execute("SELECT COUNT(*) FROM Domain").fetchone()[0]
    print(f"3did DDI1 rows after:   {after_ddi}", flush=True)
    print(f"3did Domain rows after: {after_dom}", flush=True)
    con.close()

    con = sqlite3.connect("3did.smoke.sqlite3")
    con.execute("VACUUM")
    con.close()

    # ---- Negatome: keep N random non-empty lines
    with open("${negatome}") as fh:
        lines = [ln for ln in fh if ln.strip()]
    print(f"negatome rows before: {len(lines)}", flush=True)
    random.shuffle(lines)
    sampled = lines[:min(N, len(lines))]
    with open("negatome.smoke.tsv", "w") as fh:
        fh.writelines(sampled)
    print(f"negatome rows after:  {len(sampled)}", flush=True)

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\\n")
    """
}
