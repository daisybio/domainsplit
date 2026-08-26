process INSERT_EXTERNAL_SOURCES {
    tag "insert_external_sources"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path external_ddis

    output:
    path "domainsplit.sqlite3",   emit: domainsplit_db
    path "source_conflicts.tsv",  emit: conflicts
    path "versions.yml",          emit: versions

    script:
    """
    cp --reflink=auto "${domainsplit_db_in}" domainsplit.sqlite3

    insert_external_sources.py \\
        --db domainsplit.sqlite3 \\
        --ddis ${external_ddis} \\
        --conflicts-out source_conflicts.tsv \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
