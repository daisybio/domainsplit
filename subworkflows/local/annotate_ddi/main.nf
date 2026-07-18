/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SCORE DDI - Calculate z-scores for domain-domain interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Workflow overview
-----------------

1.  BUILD_SCORING_MATRIX
        Derive empirical potential based on the positive interactions in the training set

2.  SCORE_DDI
        Calculates the z-score for all DDIs in training, testing and optimization split
        Scores are calculated for every DDI instance and then aggregated for each DDI.
        New column interaction_confirmed (true/false) is added to the ddi table based on the z-score threshold.
----------------------------------------------------------------------------*/


include { BUILD_SCORING_MATRIX } from '../../../modules/local/structural/build_scoring_matrix/main.nf'
include { SCORE_DDI            } from '../../../modules/local/structural/score_ddi/main.nf'


workflow ANNOTATE_DDI {

    take:
    split_db   // path to splitdatabase

    main:

    ch_versions = Channel.empty()

    // Step 1: Build scoring matrix based on 3DID empirical potential
    // only need training set here
    train_ch = split_db
        .filter { meta, db -> meta.split == 'train' }
        .map { meta, db -> [ [id: meta.method], db ] }

    BUILD_SCORING_MATRIX(train_ch)

    matrix_ch = BUILD_SCORING_MATRIX.out.c_ab_matrix
        .join(BUILD_SCORING_MATRIX.out.db_freq)
        .join(BUILD_SCORING_MATRIX.out.t_db)


    score_input_ch = split_db
        .map { meta, db -> [ [id: meta.method], meta, db ] }
        .combine(matrix_ch, by: 0)
        .map { method_meta, meta, db, c_ab_matrix, db_freq, t_db ->
            [ meta, db, c_ab_matrix, db_freq, t_db ]
        }

    // Step 2: Assign z-scores to all ddi complexes in the database based on the empirical potential
    SCORE_DDI(
        score_input_ch
    )

    ch_versions = ch_versions.mix(
        BUILD_SCORING_MATRIX.out.versions,
        SCORE_DDI.out.versions,
    )
        

    emit:
    scores   = SCORE_DDI.out.dbscored  // [ meta, scored dbscored ]
    versions = ch_versions
}
