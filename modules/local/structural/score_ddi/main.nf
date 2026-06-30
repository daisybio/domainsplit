process SCORE_DDI {
    tag "score_ddi"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path db, stageAs: 'input.domainsplit.sqlite3'
    path c_ab_matrix, stageAs: 'c_ab_matrix.csv'
    path db_freq,     stageAs: 'db_freq.csv'
    path t_db,       stageAs: 't_db.csv'

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:

    """
    score_ddi.py \\
        --db_in  ${db} \\
        --db_out domainsplit.sqlite3 \\
        --zscore_threshold ${params.zscore_threshold ?: 2.3} \\
        --c_ab_matrix ${c_ab_matrix} \\
        --db_freq ${db_freq} \\
        --t_db ${t_db} \\
        --versions versions.yml \\
        --process_name "${task.process}"
    """
}
