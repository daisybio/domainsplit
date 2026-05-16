process DOWNLOAD_3DID_SQLITE {
    tag "3did"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    // No container: needs awk + zcat (host) + python3 stdlib sqlite3.
    // Host has all three on the cluster compute nodes; avoids needing a
    // combo container (no quay image bundles sqlite + gawk + gzip).

    input:
    path mysql_gz_file

    output:
    path "3did.sqlite3", emit: sqlite
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    zcat ${mysql_gz_file} | mysql2sqlite - > 3did.dump.sql

    python3 - <<'PY'
import sqlite3
con = sqlite3.connect("3did.sqlite3")
with open("3did.dump.sql", "r") as fh:
    con.executescript(fh.read())
con.commit()
con.close()
PY

    rm -f 3did.dump.sql

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
        sqlite3: \$(python3 -c 'import sqlite3; print(sqlite3.sqlite_version)')
    END_VERSIONS
    """
}

process EXTRACT_PFAM_IDS {
    tag { "${sqlite_file.simpleName}" }
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    // No container: uses python3 stdlib sqlite3, available on host.

    input:
    path sqlite_file

    output:
    path "3did_pfam_ids.txt", emit: pfam_ids
    path "versions.yml", emit: versions

    script:
    """
    python3 - <<'PY'
import sqlite3
con = sqlite3.connect("${sqlite_file}")
cur = con.execute("SELECT DISTINCT SUBSTR(Pfam_id, 1, INSTR(Pfam_id, '.') - 1) FROM Domain;")
with open("3did_pfam_ids.txt", "w") as fh:
    for (pid,) in cur:
        if pid:
            fh.write(pid + "\\n")
con.close()
PY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
        sqlite3: \$(python3 -c 'import sqlite3; print(sqlite3.sqlite_version)')
    END_VERSIONS
    """
}
