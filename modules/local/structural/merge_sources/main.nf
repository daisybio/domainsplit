process MERGE_SOURCES {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    tuple val(meta), val(sources), path(dbs)

    output:
    tuple val(meta), path("dbscored_merged.sqlite3"), emit: dbscored
    path "versions.yml",                              emit: versions

    script:
    def db_args = [sources, dbs].transpose().collect { s, d -> "--db ${s}=${d}" }..join(' ')
    """
    merge_sources.py \\
        ${db_args} \\
        --db_out dbscored_merged.sqlite3 \\
        --versions versions.yml \\
        --process_name "${task.process}"
    """
}
