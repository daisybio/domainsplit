process INSERT_DOMAIN_GO_TERMS {
    tag "insert_domain_go_terms"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path pfam2go

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    """
    insert_domain_go_terms.py \\
        --db-in ${domainsplit_db_in} \\
        --pfam2go ${pfam2go} \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
