process INSERT_SINGLE_DOMAIN_PPI {
    tag "insert_single_domain_ppi"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path hippie_tsv
    path swissprot_map
    val  min_score

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    """
    cp "${domainsplit_db_in}" domainsplit.sqlite3

    insert_single_domain_ppi.py \\
        --db domainsplit.sqlite3 \\
        --hippie ${hippie_tsv} \\
        --swissprot-map ${swissprot_map} \\
        --min-score ${min_score} \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """

    stub:
    """
    touch domainsplit.sqlite3
    echo '"${task.process}":' > versions.yml
    echo '    stub: "true"' >> versions.yml
    """
}
