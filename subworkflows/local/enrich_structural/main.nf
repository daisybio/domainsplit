/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ENRICH_STRUCTURAL -- Process AF-derived structures and store them in a 
    HDF5 file, for later use in DDI scoring.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   Function is a wrapper around the ENRICH_STRUCTURAL_AF module,
   which does the actual work.
----------------------------------------------------------------------------*/


include { ENRICH_STRUCTURAL_AF } from '../../../modules/local/structural/enrich_structural_af/main.nf'

workflow ENRICH_STRUCTURAL {
    take:
        domainsplit_db_in
        af_metaddata

    main:
        ENRICH_STRUCTURAL_AF(
            domainsplit_db_in,
            af_metaddata
        )

        domainsplit_db = ENRICH_STRUCTURAL_AF.out.dbstruct
        structures = ENRICH_STRUCTURAL_AF.out.structures
        versions = ENRICH_STRUCTURAL_AF.out.versions

    emit:
        domainsplit_db
        structures
        versions
}