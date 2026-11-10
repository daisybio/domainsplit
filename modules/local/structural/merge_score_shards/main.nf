process MERGE_SCORE_SHARDS {
    tag "$method"
    label 'process_low'

    input:
    tuple val(method), path(score_shards, stageAs: 'scores_??.tsv'), path(confirmed_shards, stageAs: 'confirmed_??.tsv')

    output:
    tuple val(method), path("${method}.scores.tsv"), path("${method}.confirmed.tsv"), emit: merged

    script:
    """
    head -n 1 ${score_shards[0]} > ${method}.scores.tsv
    tail -n +2 -q ${score_shards} >> ${method}.scores.tsv

    head -n 1 ${confirmed_shards[0]} > ${method}.confirmed.tsv
    tail -n +2 -q ${confirmed_shards} >> ${method}.confirmed.tsv
    """
}