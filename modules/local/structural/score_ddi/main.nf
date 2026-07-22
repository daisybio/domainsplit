process SCORE_DDI {
    tag "${meta.id}:${meta.source}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    tuple val(meta), path(db, stageAs: 'input.dbsplit.sqlite3'), path(c_ab_matrix, stageAs: 'c_ab_matrix.csv'), path(db_freq, stageAs: 'db_freq.csv'), path(t_db, stageAs: 't_db.csv')

    output:
    tuple val(meta), path("dbscored.sqlite3"), emit: dbscored
    path "versions.yml",                       emit: versions

    script:
    """
    score_ddi.py \\
        --db_in  ${db} \\
        --db_out dbscored.sqlite3 \\
        --c_ab_matrix ${c_ab_matrix} \\
        --db_freq ${db_freq} \\
        --t_db ${t_db} \\
        --source ${meta.source} \\
        --versions versions.yml \\
        --process_name "${task.process}"
    """
}