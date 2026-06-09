process INSERT_PPIDM {
    tag "insert_ppidm"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path ppidm_tsv
    val  classes

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    """
    cp "${domainsplit_db_in}" domainsplit.sqlite3

    insert_ppidm.py \\
        --db domainsplit.sqlite3 \\
        --ppidm ${ppidm_tsv} \\
        --classes "${classes}" \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
