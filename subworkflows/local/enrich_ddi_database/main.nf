/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ENRICH_DDI_DATABASE -- annotate the DDI database with domain GO, protein
    sequences, protein GO and STRING PPI.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Each INSERT_* process opens the SQLite emitted by the previous step,
    performs one phase, commits, and emits the database forward. The chain
    keeps schema/transaction boundaries identical to the previous monolithic
    ENRICH process.

    Embeddings are not part of this chain. They are domain-level, published as
    HDF5 files by generate_domain_embeddings, and never enter the database -- so
    INSERT_DOMAIN_PROTEIN_MAPPING, which existed to write them onto the
    domain_protein_map rows INGEST_INSTANCES had already created, is gone with
    them: without the embedding columns it was a no-op duplicate of ingest.
----------------------------------------------------------------------------*/

include { INSERT_DOMAIN_GO_TERMS   } from '../../../modules/local/enrich/insert_domain_go_terms/main.nf'
include { INSERT_PROTEIN_SEQUENCES } from '../../../modules/local/enrich/insert_protein_sequences/main.nf'
include { INSERT_PROTEIN_GO_TERMS  } from '../../../modules/local/enrich/insert_protein_go_terms/main.nf'
include { INSERT_PPI               } from '../../../modules/local/enrich/insert_ppi/main.nf'

workflow ENRICH_DDI_DATABASE {
    take:
    domainsplit_db_in
    input_pfam2go
    input_uniprot_sequences
    protein_domain_map
    input_uniprot_go_terms
    input_string
    input_string_map

    main:
    db_after_domain_go = INSERT_DOMAIN_GO_TERMS(
        domainsplit_db_in,
        input_pfam2go
    ).domainsplit_db

    db_after_proteins = INSERT_PROTEIN_SEQUENCES(
        db_after_domain_go,
        input_uniprot_sequences,
        protein_domain_map
    ).domainsplit_db

    db_after_protein_go = INSERT_PROTEIN_GO_TERMS(
        db_after_proteins,
        input_uniprot_go_terms
    ).domainsplit_db

    domainsplit_db = INSERT_PPI(
        db_after_protein_go,
        input_string,
        input_string_map
    ).domainsplit_db

    ch_versions = Channel.empty().mix(
        INSERT_DOMAIN_GO_TERMS.out.versions,
        INSERT_PROTEIN_SEQUENCES.out.versions,
        INSERT_PROTEIN_GO_TERMS.out.versions,
        INSERT_PPI.out.versions,
    )

    emit:
    domainsplit_db
    versions = ch_versions
}
