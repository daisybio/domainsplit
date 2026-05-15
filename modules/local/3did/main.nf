process DOWNLOAD_3DID_SQLITE {
    conda "${moduleDir}/environment.yml"

    input:
    path mysql_gz_file

    output:
    path "3did.sqlite3"

    script:
    """
    zcat ${mysql_gz_file} |
    mysql2sqlite - |
    sqlite3 3did.sqlite3
    """
}

process EXTRACT_PFAM_IDS {
    conda "${moduleDir}/environment.yml"

    input:
    path sqlite_file

    output:
    path "3did_pfam_ids.txt"

    script:
    """
    sqlite3 ${sqlite_file} 'SELECT DISTINCT SUBSTR(Pfam_id, 1, INSTR(Pfam_id, ".") - 1) FROM Domain;' > 3did_pfam_ids.txt
    """
}
