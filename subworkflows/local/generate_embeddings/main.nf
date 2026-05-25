/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENERATE_EMBEDDINGS -- run protein-level ProtT5 and ESM (protein + domain)
    embedding generation in parallel branches.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ProtT5 path is driven by `params.download_prott5_complete`:
      true  -> single complete EBI HDF5 download (default; reliable but large).
      false -> chunked async REST queries against the UniProt embedding service.
    ESM path produces both per-residue protein embeddings and pooled domain
    embeddings against the supplied protein <-> domain map.
----------------------------------------------------------------------------*/

include { download_prott5_embeddings } from '../../../modules/local/prott5_embeddings/main.nf'
include { generate_esm_embeddings    } from '../../../modules/local/esm_embeddings/main.nf'

workflow GENERATE_EMBEDDINGS {
    take:
    protein_domain_map
    input_uniprot_sequences

    main:
    protein_ids = protein_domain_map.splitCsv(decompress: true, header: true)
        .map { it.uniprot_id }
        .toSortedList()
        .flatten()
        .distinct()

    download_prott5_embeddings(protein_ids)
    generate_esm_embeddings(input_uniprot_sequences, protein_domain_map)

    emit:
    prott5_embeddings       = download_prott5_embeddings.out.prott5_embeddings
    esm_protein_embeddings  = generate_esm_embeddings.out.protein_embeddings
    esm_domain_embeddings   = generate_esm_embeddings.out.domain_embeddings
}
