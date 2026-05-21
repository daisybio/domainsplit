/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ENRICH_DDI_DATABASE -- annotate the DDI database with domain GO, proteins
    + per-residue embeddings, protein GO, STRING PPI, and the per-domain
    alignment + embedding map.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Each INSERT_* process opens the SQLite emitted by the previous step,
    performs one phase, commits, and emits the database forward. The chain
    keeps schema/transaction boundaries identical to the previous monolithic
    ENRICH process.
----------------------------------------------------------------------------*/

include { INSERT_DOMAIN_GO_TERMS          } from '../../../modules/local/enrich/insert_domain_go_terms/main.nf'
include { INSERT_PROTEINS_WITH_EMBEDDINGS } from '../../../modules/local/enrich/insert_proteins_with_embeddings/main.nf'
include { INSERT_PROTEIN_GO_TERMS         } from '../../../modules/local/enrich/insert_protein_go_terms/main.nf'
include { INSERT_PPI                      } from '../../../modules/local/enrich/insert_ppi/main.nf'
include { INSERT_DOMAIN_PROTEIN_MAPPING   } from '../../../modules/local/enrich/insert_domain_protein_mapping/main.nf'

workflow ENRICH_DDI_DATABASE {
    take:
    domainsplit_db_in
    input_pfam2go
    input_uniprot_sequences
    protein_domain_map
    prott5_embeddings
    input_uniprot_go_terms
    input_string
    input_uniprot_id_mapping
    esm_protein_embeddings
    esm_domain_embeddings

    main:
    db_after_domain_go = INSERT_DOMAIN_GO_TERMS(
        domainsplit_db_in,
        input_pfam2go
    ).domainsplit_db

    db_after_proteins = INSERT_PROTEINS_WITH_EMBEDDINGS(
        db_after_domain_go,
        input_uniprot_sequences,
        protein_domain_map,
        prott5_embeddings,
        esm_protein_embeddings
    ).domainsplit_db

    db_after_protein_go = INSERT_PROTEIN_GO_TERMS(
        db_after_proteins,
        input_uniprot_go_terms
    ).domainsplit_db

    db_after_ppi = INSERT_PPI(
        db_after_protein_go,
        input_string,
        input_uniprot_id_mapping
    ).domainsplit_db

    domainsplit_db = INSERT_DOMAIN_PROTEIN_MAPPING(
        db_after_ppi,
        protein_domain_map,
        esm_domain_embeddings
    ).domainsplit_db.first()

    emit:
    domainsplit_db
}
