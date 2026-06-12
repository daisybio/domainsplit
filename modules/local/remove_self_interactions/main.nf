process REMOVE_SELF_INTERACTIONS {
    tag "remove_self_interactions"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    """
    #!/usr/bin/env python3
    import shutil
    import sqlite3
    import sys

    shutil.copy("${domainsplit_db_in}", "domainsplit.sqlite3")

    conn = sqlite3.connect("domainsplit.sqlite3")
    conn.execute("PRAGMA foreign_keys=ON")
    before = conn.execute("SELECT COUNT(*) FROM domain_domain_interaction").fetchone()[0]
    conn.execute(
        "DELETE FROM domain_domain_interaction WHERE domain_id_a = domain_id_b"
    )
    conn.commit()
    after = conn.execute("SELECT COUNT(*) FROM domain_domain_interaction").fetchone()[0]
    conn.close()
    print(f"[remove_self_interactions] removed {before - after} self-DDIs "
          f"({before} -> {after})", flush=True)

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\\n")
    """

    stub:
    """
    touch domainsplit.sqlite3
    echo '"${task.process}":' > versions.yml
    echo '    stub: "true"' >> versions.yml
    """
}
