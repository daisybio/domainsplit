process INSERT_3DID {
    tag "insert_3did"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path sqlite_3did

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    """
    cp --reflink=auto "${domainsplit_db_in}" domainsplit.sqlite3

    insert_3did.py \\
        --db domainsplit.sqlite3 \\
        --sqlite-3did ${sqlite_3did} \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
