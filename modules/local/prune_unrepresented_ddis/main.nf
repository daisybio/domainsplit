process PRUNE_UNREPRESENTED_DDIS {
    tag "prune_unrepresented_ddis"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "pruned_ddis.tsv",     emit: pruned
    path "versions.yml",        emit: versions

    script:
    """
    cp --reflink=auto "${domainsplit_db_in}" domainsplit.sqlite3

    prune_unrepresented_ddis.py \\
        --db domainsplit.sqlite3 \\
        --report-out pruned_ddis.tsv \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
