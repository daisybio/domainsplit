/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CURATE_DOMAINS -- enumerate the unique Pfam domains referenced by the DDI
    set, fetch their Pfam alignments, and build a protein <-> domain map.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Reads only the post-DDI domainsplit SQLite (no 3did handle). The set of
    Pfam IDs is the projection of domain rows still referenced by
    domain_domain_interaction after smoke filtering.
----------------------------------------------------------------------------*/

include { EXTRACT_UNIQUE_DOMAINS                                  } from '../../../modules/local/curate_domains/extract_unique_domains/main.nf'
include { CREATE_PROTEIN_DOMAIN_MAPPING; DOWNLOAD_PFAM_ALIGNMENT  } from '../../../modules/local/pfam/main.nf'

workflow CURATE_DOMAINS {
    take:
    domainsplit_db
    input_uniprot_id_mapping

    main:
    pfam_ids = EXTRACT_UNIQUE_DOMAINS(domainsplit_db).pfam_ids.splitText()

    pfam_files = DOWNLOAD_PFAM_ALIGNMENT(pfam_ids).alignment

    protein_domain_map = CREATE_PROTEIN_DOMAIN_MAPPING(
        input_uniprot_id_mapping,
        pfam_files.collect()
    ).mapping

    emit:
    protein_domain_map
}
