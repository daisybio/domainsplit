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
include { generate_esm_embeddings     } from '../modules/local/esm_embeddings/main.nf'
include { ENRICH_DDI_DATABASE         } from '../subworkflows/local/enrich_ddi_database/main.nf'
include { ENRICH_STRUCTURAL           } from '../subworkflows/local/enrich_structural/main.nf'
include { SPLIT_DOMAINSPLIT_DATABASE  } from '../subworkflows/local/split_domainsplit_database/main.nf'
include { ANNOTATE_DDI                   } from '../subworkflows/local/annotate_ddi/main.nf'
include { ANALYZE_DDI_BIAS            } from '../modules/local/analyze_ddi_bias/main.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DOMAINSPLIT {
main:
    ch_versions = Channel.empty()

    input_uniprot_id_mapping = file(params.url_uniprot_id_mapping)
    input_uniprot_go_terms   = file(params.url_uniprot_go_terms)
    input_uniprot_sequences  = file(params.url_uniprot_sequences)
    input_string             = file(params.url_string)
    input_pfam2go            = file(params.url_pfam2go)
    input_3did               = file(params.url_3did)

    def prott5_file = file(params.url_uniprot_prott5_embeddings)

    empty_db = INIT_DOMAINSPLIT_DB().domainsplit_db

    COLLECT_DDI_DATA(
        empty_db,
        params.url_3did,
        params.url_negatome,
        params.url_uniprot_swissprot_pfam,
        params.hippie_tsv,
        params.ppidm_tsv,
        params.negative_ppi_parquet,
    )

    domainsplit_db_ddi = COLLECT_DDI_DATA.out.domainsplit_db

    CURATE_DOMAINS(
        domainsplit_db_ddi,
        input_uniprot_id_mapping,
    )

    protein_domain_map = CURATE_DOMAINS.out.protein_domain_map

    generate_esm_embeddings(
        input_uniprot_sequences,
        protein_domain_map,
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
        generate_esm_embeddings.out.protein_embeddings,
        generate_esm_embeddings.out.domain_embeddings,
    )


    ENRICH_STRUCTURAL(
        ENRICH_DDI_DATABASE.out.domainsplit_db
    )

    ANALYZE_DDI_BIAS(
        ENRICH_DDI_DATABASE.out.domainsplit_db
    )

    SPLIT_DOMAINSPLIT_DATABASE(
        ENRICH_STRUCTURAL.out.domainsplit_db
    )

    ANNOTATE_DDI(
        SPLIT_DOMAINSPLIT_DATABASE.out.split_db
    )


    //
    // Collate and save software versions
    //
    ch_versions = ch_versions.mix(
        INIT_DOMAINSPLIT_DB.out.versions,
        COLLECT_DDI_DATA.out.versions,
        CURATE_DOMAINS.out.versions,
        generate_esm_embeddings.out.versions,
        ENRICH_DDI_DATABASE.out.versions,
        ANALYZE_DDI_BIAS.out.versions,
        SPLIT_DOMAINSPLIT_DATABASE.out.versions,
        ANNOTATE_DDI.out.versions
    )

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'pipeline_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

emit:
    domainsplit_db  = ENRICH_DDI_DATABASE.out.domainsplit_db
    split_db        = SPLIT_DOMAINSPLIT_DATABASE.out.split_db
    scored_db       = ANNOTATE_DDI.out.scores
    bias_report     = ANALYZE_DDI_BIAS.out.report_dir
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
