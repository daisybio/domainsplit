/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENERATE_EMBEDDINGS -- run protein-level ESM (protein + domain)
    embedding generation.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ESM path produces both per-residue protein embeddings and pooled domain
    embeddings against the supplied protein <-> domain map.

    ProtT5 embeddings are supplied externally via params.prott5_per_residue_h5
    and resolved in the top-level workflow (domainsplit.nf).
----------------------------------------------------------------------------*/

include { generate_esm_embeddings    } from '../../../modules/local/esm_embeddings/main.nf'

workflow GENERATE_EMBEDDINGS {
    take:
    protein_domain_map
    input_uniprot_sequences

    main:
    generate_esm_embeddings(input_uniprot_sequences, protein_domain_map)

    emit:
    esm_protein_embeddings  = generate_esm_embeddings.out.protein_embeddings
    esm_domain_embeddings   = generate_esm_embeddings.out.domain_embeddings
}
