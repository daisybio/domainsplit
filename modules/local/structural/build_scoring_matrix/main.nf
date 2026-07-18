process BUILD_SCORING_MATRIX {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    tuple val(meta), path(dbtrain, stageAs: 'input.dbtrain.sqlite3')

    output:
    tuple val(meta), path("c_ab_matrix.csv"), emit: c_ab_matrix
    tuple val(meta), path("db_freq.csv"),     emit: db_freq
    tuple val(meta), path("t_db.txt"),        emit: t_db
    path "versions.yml",                      emit: versions

    script:
    """
    build_scoring_matrix.py \\
        --db_in ${dbtrain} \\
        --c_ab_matrix c_ab_matrix.csv \\
        --db_freq db_freq.csv \\
        --t_db t_db.txt \\
        --versions versions.yml \\
        --process_name "${task.process}"
    """
}