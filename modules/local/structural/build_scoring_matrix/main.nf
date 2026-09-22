// modules/local/structural/build_scoring_matrix/main.nf
process BUILD_SCORING_MATRIX {
    tag "$method"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    tuple val(method), path(dbstruct, stageAs: 'input.dbstruct.sqlite3'), path(structures_h5)

    output:
    tuple val(method), path("${method}.c_ab_matrix.csv"), path("${method}.db_freq.csv"), path("${method}.t_db.txt"), emit: matrix
    path "versions.yml", emit: versions

    script:
    """
    build_scoring_matrix.py \\
        --db_in input.dbstruct.sqlite3 --structures_h5 ${structures_h5} --method ${method} \\
        --c_ab_matrix ${method}.c_ab_matrix.csv --db_freq ${method}.db_freq.csv --t_db ${method}.t_db.txt \\
        --versions versions.yml --process_name "${task.process}"
    """
}