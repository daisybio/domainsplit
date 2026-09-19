// subworkflows/local/annotate_ddi/main.nf
include { BUILD_SCORING_MATRIX } from '../../../modules/local/structural/build_scoring_matrix/main.nf'
include { SCORE_DDI            } from '../../../modules/local/structural/score_ddi/main.nf'
include { MERGE_SCORES         } from '../../../modules/local/structural/merge_scores/main.nf'

workflow ANNOTATE_DDI {
    take:
    domainsplit_db
    structures

    main:
    ch_versions = Channel.empty()

    database  = domainsplit_db.first()   // broadcastable: every consumer sees the same one item
    hdf_structures = structures.first()

    methods_ch = Channel.fromList(params.structure_data_methods.toString().tokenize(','))
    .combine(database)
    .combine(hdf_structures)

    BUILD_SCORING_MATRIX(methods_ch)
    ch_versions = ch_versions.mix(BUILD_SCORING_MATRIX.out.versions)

    score_in = methods_ch.join(BUILD_SCORING_MATRIX.out.matrix, by: 0)

    SCORE_DDI(score_in)
    ch_versions = ch_versions.mix(SCORE_DDI.out.versions)

    MERGE_SCORES(
        database,
        SCORE_DDI.out.scores.map { it[1] }.collect(),
        SCORE_DDI.out.confirmed.map { it[1] }.collect(),
    )
    ch_versions = ch_versions.mix(MERGE_SCORES.out.versions)

    emit:
    domainsplit_db = MERGE_SCORES.out.dbstruct
    versions       = ch_versions
}