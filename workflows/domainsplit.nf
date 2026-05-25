/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_domainsplit_pipeline'
include { INIT_DOMAINSPLIT_DB         } from '../modules/local/init_domainsplit_db/main.nf'
include { COLLECT_DDI_DATA            } from '../subworkflows/local/collect_ddi_data/main.nf'
include { CURATE_DOMAINS              } from '../subworkflows/local/curate_domains/main.nf'
include { GENERATE_EMBEDDINGS         } from '../subworkflows/local/generate_embeddings/main.nf'
include { ENRICH_DDI_DATABASE         } from '../subworkflows/local/enrich_ddi_database/main.nf'
include { SPLIT_DOMAINSPLIT_DATABASE  } from '../subworkflows/local/split_domainsplit_database/main.nf'
include { ANALYZE_DDI_BIAS            } from '../modules/local/analyze_ddi_bias/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DOMAINSPLIT {
main:
    input_uniprot_id_mapping = file(params.url_uniprot_id_mapping)
    input_uniprot_embeddings = file(params.url_uniprot_embeddings)
    input_uniprot_go_terms   = file(params.url_uniprot_go_terms)
    input_uniprot_sequences  = file(params.url_uniprot_sequences)
    input_string             = file(params.url_string)
    input_pfam2go            = file(params.url_pfam2go)

    def prott5_file = []
    if (params.prott5_per_residue_h5) {
        def f = file(params.prott5_per_residue_h5)
        if (f.exists()) {
            prott5_file = f
        } else {
            log.warn "ProtT5 HDF5 not found at '${params.prott5_per_residue_h5}' — skipping ProtT5 embeddings"
        }
    } else {
        log.warn "params.prott5_per_residue_h5 not set — skipping ProtT5 embeddings"
    }

    empty_db = INIT_DOMAINSPLIT_DB().domainsplit_db

    def pfam_stockholm_file = params.pfam_stockholm
        ? file(params.pfam_stockholm, checkIfExists: true)
        : []

    COLLECT_DDI_DATA(
        empty_db,
        params.url_3did,
        params.url_negatome,
        input_uniprot_id_mapping,
        pfam_stockholm_file,
    )

    domainsplit_db_ddi = COLLECT_DDI_DATA.out.domainsplit_db

    CURATE_DOMAINS(
        domainsplit_db_ddi,
        input_uniprot_id_mapping,
    )

    protein_domain_map = CURATE_DOMAINS.out.protein_domain_map

    GENERATE_EMBEDDINGS(
        protein_domain_map,
        input_uniprot_sequences,
    )

    ENRICH_DDI_DATABASE(
        domainsplit_db_ddi,
        input_pfam2go,
        input_uniprot_sequences,
        protein_domain_map,
        prott5_file,
        input_uniprot_go_terms,
        input_string,
        input_uniprot_id_mapping,
        GENERATE_EMBEDDINGS.out.esm_protein_embeddings,
        GENERATE_EMBEDDINGS.out.esm_domain_embeddings,
    )

    ANALYZE_DDI_BIAS(
        ENRICH_DDI_DATABASE.out.domainsplit_db
    )

    SPLIT_DOMAINSPLIT_DATABASE(
        ENRICH_DDI_DATABASE.out.domainsplit_db
    )

emit:
    domainsplit_db  = ENRICH_DDI_DATABASE.out.domainsplit_db
    split_db        = SPLIT_DOMAINSPLIT_DATABASE.out.split_db
    bias_report     = ANALYZE_DDI_BIAS.out.report_dir
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
