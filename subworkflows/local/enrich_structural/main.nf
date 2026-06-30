/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ENRICH_STRUCTURAL -- annotate the DDI database with structural information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Workflow overview
-----------------
1.  PREDICT_COMPLEX_AF / PREDICT_COMPLEX_RF
        Predict co-complex structures for all PPIs with AF3 and/or RF2.
        Each predictor slices the relevant domains from the predicted
        structure using the protein–domain mapping (res_start/res_end only;
        chain is derived from the prediction order and recorded in
        `complex_chain_map`), and stores the gzip-compressed PDB bytes in
        `domain_structure` with source='AF3' or source='RF2'.


2. DOMAIN_SLICE
        Slices every DDI in every PPI context out of the predicted co-complex
        structures and stores them in `domain_structure`.               

3.  BUILD_SCORING_MATRIX
        Derives the database-wide empirical potential (C_ab matrix,
        surface frequencies, T_DB) from the 3did structural instances 

4.  SCORE_DDI
        Calculates the z-score for every DDI in the database based on the empirical
        potential in all protein contexts where the DDI is observed.  The z-score is stored in
        domain_structure.z_score

----------------------------------------------------------------------------*/


include { BUILD_SCORING_MATRIX } from '../../../modules/local/structural/build_scoring_matrix/main.nf'
include { PREDICT_COMPLEX_RF   } from '../../../modules/local/structural/predict_complex_rf/main.nf'
include { PREDICT_COMPLEX_AF   } from '../../../modules/local/structural/predict_complex_af/main.nf'
include { DOMAIN_SLICE         } from '../../../modules/local/structural/domain_slicing/main.nf'
include { SCORE_DDI            } from '../../../modules/local/structural/score_ddi/main.nf'

workflow ENRICH_STRUCTURAL {

    take:
    domainsplit_db_in   // channel: path to SQLite DB (with ddi_complex_coords populated)

    main:

    // Step1: Predict co-complex structures for all PPIs with AF3 and/or RF2

    af_db_ch = Channel.empty()
    rf_db_ch = Channel.empty()

    if (params.alphafold) {
        PREDICT_COMPLEX_AF(domainsplit_db_in)
        af_db_ch = PREDICT_COMPLEX_AF.out.domainsplit_db
        af_ppi_db = PREDICT_COMPLEX_AF.out.pdb_files
    }

    if (params.rosettafold) {
        PREDICT_COMPLEX_RF(domainsplit_db_in)
        rf_db_ch = PREDICT_COMPLEX_RF.out.domainsplit_db
        rf_ppi_db = PREDICT_COMPLEX_RF.out.pdb_files
    }


    // Collapse to a single DB channel; fall back to the input DB if neither
    // predictor ran (so downstream steps still receive a valid database).
    db_after_prediction = af_db_ch
        .mix(rf_db_ch)
        .ifEmpty(domainsplit_db_in)
        .first()

    ppi_files = af_ppi_db
        .mix(rf_ppi_db)
        .ifEmpty(Channel.empty())
        .first()

    // Step 2: Slice domains and store them
    DOMAIN_SLICE(db_after_prediction, ppi_files)

    // Step 3: Build scoring matrix based on 3DID empirical potential
    BUILD_SCORING_MATRIX(DOMAIN_SLICE.out.domainsplit_db)

    // Step 4: Assign z-scores to all ddi complexes in the database based on the empirical potential
    SCORE_DDI(
        DOMAIN_SLICE.out.domainsplit_db,
        BUILD_SCORING_MATRIX.out.c_ab_matrix,
        BUILD_SCORING_MATRIX.out.db_freq,
        BUILD_SCORING_MATRIX.out.t_db
    )
        

    emit:
    domainsplit_db = SCORE_DDI.out.domainsplit_db
}
