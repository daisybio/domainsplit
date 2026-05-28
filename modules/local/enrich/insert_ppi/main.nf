process INSERT_PPI {
    tag "insert_ppi"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path string_ppi
    path uniprot_id_mapping

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    """
    insert_ppi.py \\
        --db-in ${domainsplit_db_in} \\
        --string-ppi ${string_ppi} \\
        --uniprot-id-mapping ${uniprot_id_mapping} \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
