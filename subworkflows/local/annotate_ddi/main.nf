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

3.  MERGE_SOURCES
        Merges the scored databases from different sources (AF3, RF2) into a single database.
        The final database contains all the scored DDIs with their z-scores and interaction_confirmed status for both sources.
----------------------------------------------------------------------------*/


include { BUILD_SCORING_MATRIX } from '../../../modules/local/structural/build_scoring_matrix/main.nf'
include { SCORE_DDI            } from '../../../modules/local/structural/score_ddi/main.nf'
include { MERGE_SOURCES        } from '../../../modules/local/structural/merge_sources/main.nf'


workflow ANNOTATE_DDI {

    take:
    split_db   // path to splitdatabase

    main:

    //NOTE: The scoring sources are hardcoded here for now, but can be made configurable in the future.
    sources_ch  = Channel.fromList(params.predictions ?: ['AF3', 'RF'])  // 
    ch_versions = Channel.empty()

    // Step1: Extract training set of each splot for calculating the empirical potential and the scoring matrix
    train_ch = split_db
        .filter { meta, db -> meta.split == 'train' }
        .combine(sources_ch)
        .map { meta, db, source ->
            [ [ id: "${meta.method}_${source}", method: meta.method, source: source ], db ]
        }

    // Step2: Build scoring matrix, i.e. the background distribution of domain-domain interactions based on the training set
    BUILD_SCORING_MATRIX(train_ch)

    // Combine outputs
    matrix_ch = BUILD_SCORING_MATRIX.out.c_ab_matrix
        .join(BUILD_SCORING_MATRIX.out.db_freq)
        .join(BUILD_SCORING_MATRIX.out.t_db)


    // Step3: Prepare input for scoring the DDIs in the split databases
    score_input_ch = split_db
        .combine(sources_ch)
        .map { meta, db, source ->
            def join_key = [ id: "${meta.method}_${source}", method: meta.method, source: source ]
            [ join_key, meta, db ]
        }
        .combine(matrix_ch, by: 0)
        .map { join_key, meta, db, c_ab_matrix, db_freq, t_db ->
            def score_meta = meta + [ source: join_key.source ]
            [ score_meta, db, c_ab_matrix, db_freq, t_db ]
        }

    // Step4: Assign z-scores to all ddi complexes in the database based on the empirical potential
    SCORE_DDI(
        score_input_ch
    )

    // Step5: Merge the scored databases from different sources into a single database
    per_split_ch = SCORE_DDI.out.dbscored
        .map { meta, db ->
            [ meta.subMap(['id', 'split', 'method']), meta.source, db ] 
            }
        .groupTuple(by: 0)

    single_ch = per_split_ch.filter { split_meta, sources, dbs -> sources.size() == 1 }
    multi_ch  = per_split_ch.filter { split_meta, sources, dbs -> sources.size() > 1 }

    MERGE_SOURCES(multi_ch)

    scored_db_ch = single_ch
        .map { split_meta, sources, dbs -> [ split_meta, dbs[0] ] }
        .mix(MERGE_SOURCES.out.dbscored)

    

    ch_versions = ch_versions.mix(
        BUILD_SCORING_MATRIX.out.versions,
        SCORE_DDI.out.versions,
        MERGE_SOURCES.out.versions,
    )
        

    emit:

    scored_db = scored_db_ch
    versions = ch_versions
}
