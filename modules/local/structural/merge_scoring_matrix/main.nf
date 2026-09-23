// modules/local/structural/merge_scoring_matrix/main.nf
//
// A true reduction (sum), unlike MERGE_SCORE_SHARDS -- can't just
// concatenate shard c_ab_counts.csv / surface_counts.csv, the counts need
// to be summed entrywise and db_freq.csv only normalized once, at the end.
process MERGE_SCORING_MATRIX {
    tag "$method"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    tuple val(method), path(c_ab_shards, stageAs: 'c_ab_??.csv'), path(surface_shards, stageAs: 'surface_??.csv')

    output:
    tuple val(method), path("${method}.c_ab_matrix.csv"), path("${method}.db_freq.csv"), path("${method}.t_db.txt"), emit: matrix

    script:
    """
    merge_scoring_matrix.py \\
        --method ${method} \\
        --c_ab_counts ${c_ab_shards} \\
        --surface_counts ${surface_shards} \\
        --c_ab_matrix_out ${method}.c_ab_matrix.csv \\
        --db_freq_out ${method}.db_freq.csv \\
        --t_db_out ${method}.t_db.txt
    """
}