// subworkflows/local/annotate_ddi/main.nf
include { BUILD_SCORING_MATRIX } from '../../../modules/local/structural/build_scoring_matrix/main.nf'
include { MERGE_SCORING_MATRIX } from '../../../modules/local/structural/merge_scoring_matrix/main.nf'
include { SCORE_DDI            } from '../../../modules/local/structural/score_ddi/main.nf'
include { MERGE_SCORE_SHARDS   } from '../../../modules/local/structural/merge_score_shards/main.nf'
include { MERGE_SCORES         } from '../../../modules/local/structural/merge_scores/main.nf'

// Number of parallel shards per method for each phase. Kept separate since
// build_scoring_matrix only sees the (smaller) train split while score_ddi
// scores every DDI the method places in any split -- they don't need the
// same degree of parallelism. Override in nextflow.config / -params-file;
// these are just starting points to tune against your trace.
params.build_scoring_matrix_shards = 2
params.score_ddi_shards            = 4

workflow ANNOTATE_DDI {
    take:
    domainsplit_db
    structures

    main:
    ch_versions = Channel.empty()

    database       = domainsplit_db.first()   // broadcastable: every consumer sees the same one item
    hdf_structures = structures.first()

    methods_ch = Channel.fromList(params.structure_data_methods.toString().tokenize(','))
        .combine(database)
        .combine(hdf_structures)
    // methods_ch: (method, database, hdf_structures)

    // ---- Phase 1: build each method's scoring matrix, sharded --------
    build_shard_in = methods_ch
        .combine(Channel.of(0..<params.build_scoring_matrix_shards))
        .map { method, db, h5, shard_id ->
            tuple(method, shard_id, params.build_scoring_matrix_shards, db, h5)
        }

    BUILD_SCORING_MATRIX(build_shard_in)
    ch_versions = ch_versions.mix(BUILD_SCORING_MATRIX.out.versions)

    // Group each method's shard outputs back together: groupTuple() on
    // (method, c_ab_counts, surface_counts) yields (method, [c_ab...], [surface...])
    matrix_shards_grouped = BUILD_SCORING_MATRIX.out.shard.groupTuple(by: 0)

    MERGE_SCORING_MATRIX(matrix_shards_grouped)
    // MERGE_SCORING_MATRIX.out.matrix: (method, c_ab_matrix, db_freq, t_db)

    // ---- Phase 2: score every DDI under each method, sharded ---------
    score_shard_in = methods_ch
        .combine(Channel.of(0..<params.score_ddi_shards))
        .map { method, db, h5, shard_id ->
            tuple(method, shard_id, params.score_ddi_shards, db, h5)
        }
        .combine(MERGE_SCORING_MATRIX.out.matrix, by: 0)
    // score_shard_in: (method, shard_id, n_shards, database, hdf_structures, c_ab_matrix, db_freq, t_db)
    // -- already matches SCORE_DDI's expected input tuple order

    SCORE_DDI(score_shard_in)
    ch_versions = ch_versions.mix(SCORE_DDI.out.versions)

    // Plain concatenation-safe grouping: score_ddi.py shards by ddi_id, so
    // no ddi_id repeats across a method's shards.
    score_shards_grouped = SCORE_DDI.out.shard.groupTuple(by: 0)

    MERGE_SCORE_SHARDS(score_shards_grouped)
    // MERGE_SCORE_SHARDS.out.merged: (method, scores, confirmed)

    MERGE_SCORES(
        database,
        MERGE_SCORE_SHARDS.out.merged.map { it[1] }.collect(),
        MERGE_SCORE_SHARDS.out.merged.map { it[2] }.collect(),
    )
    ch_versions = ch_versions.mix(MERGE_SCORES.out.versions)

    emit:
    domainsplit_db = MERGE_SCORES.out.dbstruct
    versions       = ch_versions
}
