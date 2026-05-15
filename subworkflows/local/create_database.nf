include { DOWNLOAD_3DID_SQLITE; EXTRACT_PFAM_IDS } from '../modules/local/3did/main.nf'
include { CREATE_PROTEIN_DOMAIN_MAPPING; DOWNLOAD_PFAM_ALIGNMENT } from '../modules/local/pfam/main.nf'
include { COBINET } from '../modules/local/cobinet/main.nf'
include { download_prott5_embeddings } from '../modules/local/prott5_embeddings/main.nf'
include { generate_esm_embeddings } from '../modules/local/esm_embeddings/main.nf'

workflow CREATE_COBINET_DATABASE {
    take:
    input_3did
    input_uniprot_id_mapping
    input_uniprot_embeddings
    input_uniprot_go_terms
    input_uniprot_sequences
    input_negatome
    input_string
    input_pfam2go

    main:
    sqlite_3did = DOWNLOAD_3DID_SQLITE(input_3did)
    pfam_ids = EXTRACT_PFAM_IDS(sqlite_3did).splitText()

    pfam_files = DOWNLOAD_PFAM_ALIGNMENT(pfam_ids)

    protein_domain_map = CREATE_PROTEIN_DOMAIN_MAPPING(input_uniprot_id_mapping, pfam_files.collect())

    protein_ids = protein_domain_map.splitCsv(decompress: true, header:true)
        .map { it.uniprot_id }
        .toSortedList()
        .flatten()
        .distinct()

    prott5_embeddings = download_prott5_embeddings(protein_ids)

    (esm_protein_embeddings, esm_domain_embeddings) = generate_esm_embeddings(
        input_uniprot_sequences,
        protein_domain_map
    )

    cobinet_db_ch = COBINET(
        ddi_3did=sqlite_3did,
        ddi_negatome=input_negatome,
        pfam2go=input_pfam2go,
        uniprot_database=input_uniprot_sequences,
        protein_domain_map=protein_domain_map,
        prott5_embeddings=prott5_embeddings,
        uniprot_go_terms=input_uniprot_go_terms,
        string_ppi=input_string,
        uniprot_id_mapping=input_uniprot_id_mapping,
        esm_protein_embeddings=esm_protein_embeddings,
        esm_domain_embeddings=esm_domain_embeddings
    ).first()

    emit:
    cobinet_db_ch
}