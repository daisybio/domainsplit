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

# Python sqlite3.executescript() forwards the entire string to sqlite3_exec(),
# bounded by SQLITE_MAX_SQL_LENGTH (default 1 MB). The 3did dump exceeds this,
# so we stream-execute one statement at a time.
#
# mysql2sqlite converts MySQL \' → '' and \\ → \ before output, so SQLite SQL
# contains no backslash escapes. We only track '' (doubled single-quote) and
# "" (doubled double-quote) — both are consumed by toggling the open-quote flag
# twice in sequence, which correctly cancels out.

con = sqlite3.connect("3did.sqlite3")

buf = []
in_squote = False
in_dquote = False

def flush():
    stmt = "".join(buf).strip()
    if not stmt:
        return
    con.executescript(stmt)

with open("3did.dump.sql", "r", encoding="utf-8", errors="replace") as fh:
    for line in fh:
        buf.append(line)
        for c in line:
            if c == "'" and not in_dquote:
                in_squote = not in_squote
            elif c == '"' and not in_squote:
                in_dquote = not in_dquote
        if not in_squote and not in_dquote and line.rstrip().endswith(";"):
            flush()
            buf = []

flush()
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
