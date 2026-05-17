process DOWNLOAD_3DID_SQLITE {
    tag "3did"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path mysql_gz_file

    output:
    path "3did.sqlite3", emit: sqlite
    path "versions.yml", emit: versions

    script:
    """
    set -euo pipefail

    zcat ${mysql_gz_file} | mysql2sqlite - > 3did.dump.sql

    3did_import.py 3did.dump.sql 3did.sqlite3

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
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path sqlite_file

    output:
    path "3did_pfam_ids.txt", emit: pfam_ids
    path "versions.yml", emit: versions

    script:
    """
    3did_extract_pfam_ids.py ${sqlite_file} 3did_pfam_ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | awk '{print \$2}')
        sqlite3: \$(python3 -c 'import sqlite3; print(sqlite3.sqlite_version)')
    END_VERSIONS
    """
}
