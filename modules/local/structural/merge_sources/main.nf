process MERGE_SOURCES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    tuple val(meta), path(db_af3, stageAs: 'af3.sqlite3'), path(db_rf, stageAs: 'rf.sqlite3')

    output:
    tuple val(meta), path("dbscored_merged.sqlite3"), emit: dbscored
    path "versions.yml",                              emit: versions

    script:
    """
    merge_sources.py \\
        --db_af3 ${db_af3} \\
        --db_rf  ${db_rf} \\
        --db_out dbscored_merged.sqlite3 \\
        --versions versions.yml \\
        --process_name "${task.process}"
    """
}
