/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_domainsplit_pipeline'
include { CREATE_COBINET_DATABASE } from './subworkflows/create_database.nf'
include { SPLIT_COBINET_DATABASE } from './subworkflows/split_database.nf'

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

    cobinet_db_ch = CREATE_COBINET_DATABASE(
        input_3did,
        input_uniprot_id_mapping,
        input_uniprot_embeddings,
        input_uniprot_go_terms,
        input_uniprot_sequences,
        input_negatome,
        input_string,
        input_pfam2go,
    )//*/

    //cobinet_db_ch = file("results/cobinet.sqlite3")

    split_db_ch = SPLIT_COBINET_DATABASE(
        cobinet_db_ch
    )

    publish:
    cobinet_db = cobinet_db_ch
    split_db = split_db_ch
}

output {
    cobinet_db  {
    }
    // put the split databases in a separate output folders
    split_db  {
        path {
            it[1] >> "split_databases/${it[0].method}/${it[0].split}.sqlite3"
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
