process DOWNLOAD_3DID_SQLITE {
    tag "3did"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    path mysql_gz_file

    output:
    path "3did.sqlite3", emit: sqlite
    path "versions.yml", emit: versions

    script:
    """
    zcat ${mysql_gz_file} |
    mysql2sqlite - |
    sqlite3 3did.sqlite3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqlite3: \$(sqlite3 --version | awk '{print \$1}')
    END_VERSIONS
    """
}

process EXTRACT_PFAM_IDS {
    tag { "${sqlite_file.simpleName}" }
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    path sqlite_file

    output:
    path "3did_pfam_ids.txt", emit: pfam_ids
    path "versions.yml", emit: versions

    script:
    """
    sqlite3 ${sqlite_file} 'SELECT DISTINCT SUBSTR(Pfam_id, 1, INSTR(Pfam_id, ".") - 1) FROM Domain;' > 3did_pfam_ids.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqlite3: \$(sqlite3 --version | awk '{print \$1}')
    END_VERSIONS
    """
}
