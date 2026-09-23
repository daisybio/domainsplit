// modules/local/structural/merge_score_shards/main.nf
process MERGE_SCORE_SHARDS {
    tag "$method"
    label 'process_low'

    input:
    tuple val(method), path(score_shards, stageAs: 'scores_??.tsv'), path(confirmed_shards, stageAs: 'confirmed_??.tsv')

    output:
    tuple val(method), path("${method}.scores.tsv"), path("${method}.confirmed.tsv"), emit: merged

    // Plain concatenation is correct here: SCORE_DDI shards by ddi_id (see
    // score_ddi.py), so ddi_ids never repeat across shards -- no shard's
    // confirmed.tsv needs re-aggregating against another shard's rows.

    script:
    """
    head -n 1 ${score_shards[0]} > ${method}.scores.tsv
    tail -n +2 -q ${score_shards} >> ${method}.scores.tsv

    head -n 1 ${confirmed_shards[0]} > ${method}.confirmed.tsv
    tail -n +2 -q ${confirmed_shards} >> ${method}.confirmed.tsv
    """
}