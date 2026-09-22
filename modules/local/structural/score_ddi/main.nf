// modules/local/structural/score_ddi/main.nf
process SCORE_DDI {
    tag "$method"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    tuple val(method), path(dbstruct, stageAs: 'input.dbstruct.sqlite3'), path(structures_h5), path(c_ab_matrix), path(db_freq), path(t_db)

    output:
    tuple val(method), path("${method}.scores.tsv"),    emit: scores
    tuple val(method), path("${method}.confirmed.tsv"), emit: confirmed
    path "versions.yml", emit: versions

    script:
    """
    score_ddi.py \\
        --db_in input.dbstruct.sqlite3 --structures_h5 ${structures_h5} --method ${method} \\
        --c_ab_matrix ${c_ab_matrix} --db_freq ${db_freq} --t_db ${t_db} \\
        --scores_out ${method}.scores.tsv --confirmed_out ${method}.confirmed.tsv \\
        --versions versions.yml --process_name "${task.process}"
    """
}