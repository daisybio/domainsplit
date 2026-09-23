// modules/local/structural/build_scoring_matrix/main.nf
process BUILD_SCORING_MATRIX {
    tag "$method-$shard_id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    tuple val(method), val(shard_id), val(n_shards),
          path(dbstruct, stageAs: 'input.dbstruct.sqlite3'), path(structures_h5)

    output:
    tuple val(method), path("${method}.${shard_id}.c_ab_counts.csv"),
                        path("${method}.${shard_id}.surface_counts.csv"), emit: shard
    path "versions.yml", emit: versions

    script:
    """
    build_scoring_matrix.py \\
        --db_in input.dbstruct.sqlite3 --structures_h5 ${structures_h5} --method ${method} \\
        --c_ab_counts_out ${method}.${shard_id}.c_ab_counts.csv \\
        --surface_counts_out ${method}.${shard_id}.surface_counts.csv \\
        --shard_id ${shard_id} --n_shards ${n_shards} \\
        --versions versions.yml --process_name "${task.process}"
    """
}