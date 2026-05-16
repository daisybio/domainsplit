/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CREATE_COBINET_DATABASE -- assemble the unified CoBiNet SQLite database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Downloads and normalises 3did, Pfam, ProtT5, and ESM inputs, then
    assembles them into a single cobinet.sqlite3 file via the COBINET module.
----------------------------------------------------------------------------*/

include { DOWNLOAD_3DID_SQLITE; EXTRACT_PFAM_IDS                  } from '../../../modules/local/3did/main.nf'
include { CREATE_PROTEIN_DOMAIN_MAPPING; DOWNLOAD_PFAM_ALIGNMENT  } from '../../../modules/local/pfam/main.nf'
include { COBINET                                                 } from '../../../modules/local/cobinet/main.nf'
include { download_prott5_embeddings                              } from '../../../modules/local/prott5_embeddings/main.nf'
include { generate_esm_embeddings                                 } from '../../../modules/local/esm_embeddings/main.nf'

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
    sqlite_3did = DOWNLOAD_3DID_SQLITE(input_3did).sqlite
    pfam_ids = EXTRACT_PFAM_IDS(sqlite_3did).pfam_ids.splitText()

    pfam_files = DOWNLOAD_PFAM_ALIGNMENT(pfam_ids).alignment

    protein_domain_map = CREATE_PROTEIN_DOMAIN_MAPPING(input_uniprot_id_mapping, pfam_files.collect()).mapping

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
        sqlite_3did,
        input_negatome,
        input_pfam2go,
        input_uniprot_sequences,
        protein_domain_map,
        prott5_embeddings,
        input_uniprot_go_terms,
        input_string,
        input_uniprot_id_mapping,
        esm_protein_embeddings,
        esm_domain_embeddings
    ).cobinet_db.first()

    emit:
    cobinet_db = cobinet_db_ch
}
