/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_domainsplit_pipeline'
include { CREATE_COBINET_DATABASE } from '../subworkflows/local/create_cobinet_database/main.nf'
include { SPLIT_COBINET_DATABASE  } from '../subworkflows/local/split_cobinet_database/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DOMAINSPLIT {
main:
    input_3did               = file(params.url_3did)
    input_uniprot_id_mapping = file(params.url_uniprot_id_mapping)
    input_uniprot_embeddings = file(params.url_uniprot_embeddings)
    input_uniprot_go_terms   = file(params.url_uniprot_go_terms)
    input_uniprot_sequences  = file(params.url_uniprot_sequences)
    input_negatome           = file(params.url_negatome)
    input_string             = file(params.url_string)
    input_pfam2go            = file(params.url_pfam2go)

    CREATE_COBINET_DATABASE(
        input_3did,
        input_uniprot_id_mapping,
        input_uniprot_embeddings,
        input_uniprot_go_terms,
        input_uniprot_sequences,
        input_negatome,
        input_string,
        input_pfam2go,
    )

    SPLIT_COBINET_DATABASE(
        CREATE_COBINET_DATABASE.out.cobinet_db
    )

emit:
    cobinet_db = CREATE_COBINET_DATABASE.out.cobinet_db
    split_db   = SPLIT_COBINET_DATABASE.out.split_db
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
