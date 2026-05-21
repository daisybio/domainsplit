process INSERT_DDIS {
    tag "insert_ddis"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path cobinet_db_in, stageAs: 'input.cobinet.sqlite3'
    path sqlite_3did
    path negatome_txt

    output:
    path "cobinet.sqlite3", emit: cobinet_db
    path "versions.yml",    emit: versions

    script:
    """
    #!/usr/bin/env python3
    import shutil
    import sqlite3
    import sys

    shutil.copy("${cobinet_db_in}", "cobinet.sqlite3")

    conn_3did    = sqlite3.connect("${sqlite_3did}")
    conn_cobinet = sqlite3.connect("cobinet.sqlite3")
    conn_cobinet.execute("PRAGMA foreign_keys=ON")
    conn_cobinet.execute("PRAGMA journal_mode=OFF")
    conn_cobinet.execute("PRAGMA synchronous=OFF")

    def iter_negatome_pairs(path):
        with open(path) as f:
            for line in f:
                tokens = line.split()
                if len(tokens) < 2:
                    continue
                yield tokens[0], tokens[1]

    def negatome_pfam_ids(path):
        ids = set()
        for a, b in iter_negatome_pairs(path):
            ids.add(a)
            ids.add(b)
        return ids

    # ---- domain rows: 3did Domain x domain_length + negatome pfam ids
    print("Inserting domain information", flush=True)
    cursor = conn_3did.execute(
        "SELECT Name, Pfam_id, profile_length "
        "FROM Domain, domain_length "
        "WHERE domain_length.domain = Domain.Name"
    )
    domain_rows_3did = ((name, pfam_id.split(".")[0]) for (name, pfam_id, _length) in cursor)
    domain_rows_negatome = ((None, pfam_id) for pfam_id in negatome_pfam_ids("${negatome_txt}"))

    conn_cobinet.executemany(
        "INSERT OR IGNORE INTO domain(name, pfam_id) VALUES (?, ?);",
        list(domain_rows_3did) + list(domain_rows_negatome),
    )
    cursor.close()
    conn_cobinet.commit()

    # ---- positive DDIs from 3did
    print("Inserting positive DDIs from 3did", flush=True)
    cursor = conn_3did.execute(
        "SELECT d1.Pfam_id, d2.Pfam_id "
        "FROM DDI1, Domain AS d1, Domain AS d2 "
        "WHERE DDI1.domain1 = d1.Name AND DDI1.domain2 = d2.Name;"
    )
    pos_iter = ((id_1.split(".")[0], id_2.split(".")[0]) for (id_1, id_2) in cursor)
    conn_cobinet.executemany(
        '''INSERT OR IGNORE INTO domain_domain_interaction(domain_id_a, domain_id_b, negative, source)
           SELECT d1.id, d2.id, FALSE, '3did'
           FROM domain AS d1, domain AS d2
           WHERE d1.pfam_id = ? AND d2.pfam_id = ?;''',
        pos_iter,
    )
    cursor.close()
    conn_cobinet.commit()

    # ---- negative DDIs from negatome
    print("Inserting negative DDIs from negatome", flush=True)
    conn_cobinet.executemany(
        '''INSERT OR IGNORE INTO domain_domain_interaction(domain_id_a, domain_id_b, negative, source)
           SELECT d1.id, d2.id, TRUE, 'negatome'
           FROM domain AS d1, domain AS d2
           WHERE d1.pfam_id = ? AND d2.pfam_id = ?;''',
        iter_negatome_pairs("${negatome_txt}"),
    )
    conn_cobinet.commit()
    conn_cobinet.close()
    conn_3did.close()

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\\n")
    """
}
