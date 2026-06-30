process BUILD_SCORING_MATRIX {
    tag "build_scoring_matrix"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0" 

    input:
        path db

    output:
        // c_ab_matrix: CSV with columns chain_id_a, chain_id_b, c_ab
        path "c_ab_matrix.csv", emit: c_ab_matrix
        // db_freq: CSV with columns chain_id, frequency
        path "db_freq.csv", emit: db_freq
        // t_db: CSV with columns chain_id, t_db
        path "t_db.txt", emit: t_db

    script:
        """
        build_scoring_matrix.py \\
            --db_in ${db} \\
            --c_ab_matrix c_ab_matrix.csv \\
            --db_freq db_freq.csv \\
            --t_db t_db.txt \\
            --versions versions.yml \\
            --process_name "${task.process}"
        """
}