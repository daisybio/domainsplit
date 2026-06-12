process INSERT_NEGATOME {
    tag "insert_negatome"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path negatome_txt

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    """
    cp "${domainsplit_db_in}" domainsplit.sqlite3

    insert_negatome.py \\
        --db domainsplit.sqlite3 \\
        --negatome ${negatome_txt} \\
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
