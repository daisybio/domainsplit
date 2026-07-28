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

----------------------------------------------------------------------------*/


include { PREDICT_COMPLEX_RF   } from '../../../modules/local/structural/predict_complex_rf/main.nf'
include { PREDICT_COMPLEX_AF   } from '../../../modules/local/structural/predict_complex_af/main.nf'
include { DOMAIN_SLICE         } from '../../../modules/local/structural/domain_slicing/main.nf'


workflow ENRICH_STRUCTURAL {

    take:
    domainsplit_db_in   // channel: path to SQLite DB (with ddi_complex_coords populated)
    meta_ppis

    main:

    // Step1: Predict co-complex structures for all PPIs with AF3 and/or RF2
    af_ppi_db = Channel.empty()
    rf_ppi_db = Channel.empty()

    PREDICT_COMPLEX_AF(domainsplit_db_in, meta_ppis)
    af_ppi_db = PREDICT_COMPLEX_AF.out.pdb_files

    PREDICT_COMPLEX_RF(domainsplit_db_in, meta_ppis)
    rf_ppi_db = PREDICT_COMPLEX_RF.out.pdb_files

    // Step 2: Slice domains and store them
    DOMAIN_SLICE(domainsplit_db_in, af_ppi_db, rf_ppi_db)

        

    emit:
    domainsplit_db = DOMAIN_SLICE.out.dbstruct
}
