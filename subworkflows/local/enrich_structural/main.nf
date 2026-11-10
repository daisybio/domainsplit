/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ENRICH_STRUCTURAL -- slice every instance pair of the given methods out of
    its AF3 complex into structures.h5 and set ddi_split_membership.has_structure.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ENRICH_STRUCTURAL_AF } from '../../../modules/local/structural/enrich_structural_af/main.nf'

workflow ENRICH_STRUCTURAL {
    take:
    domainsplit_db_in // the enriched master
    af_metadata       // CSV uniprot_id_a,uniprot_id_b,path,source
    methods           // comma-separated methods whose instance pairs are sliced

    main:
    ENRICH_STRUCTURAL_AF(domainsplit_db_in, af_metadata, methods)

    emit:
    domainsplit_db = ENRICH_STRUCTURAL_AF.out.dbstruct
    structures     = ENRICH_STRUCTURAL_AF.out.structures
    versions       = ENRICH_STRUCTURAL_AF.out.versions
}
