/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT_DDI_DATA -- download and parse every DDI source, write into the
    pre-initialised Domainsplit SQLite. Downstream code consumes only the
    database; the 3did SQLite stays internal to this subworkflow.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Add a new DDI source by:
      1. Adding its download module (network fetch + format normalisation).
      2. Calling it here and routing the parsed output into INSERT_DDIS (or a
         per-source INSERT_* module if its parsing differs).
      3. Tagging its rows with a unique source string in domain_domain_interaction.
----------------------------------------------------------------------------*/

include { DOWNLOAD_3DID_SQLITE      } from '../../../modules/local/3did/main.nf'
include { DOWNLOAD_NEGATOME         } from '../../../modules/local/negatome/main.nf'
include { INSERT_DDIS               } from '../../../modules/local/insert_ddis/main.nf'
include { INSERT_PPI_NEGATIVE_DDIS  } from '../../../modules/local/insert_ppi_negative_ddis/main.nf'
include { SMOKE_FILTER              } from '../../../modules/local/smoke_filter/main.nf'

workflow COLLECT_DDI_DATA {
    take:
    domainsplit_db_in
    url_3did
    url_negatome
    uniprot_id_mapping

    main:
    file_3did     = file(url_3did)
    sqlite_3did   = DOWNLOAD_3DID_SQLITE(file_3did).sqlite
    negatome_file = DOWNLOAD_NEGATOME(url_negatome).negatome

    domainsplit_db = INSERT_DDIS(domainsplit_db_in, sqlite_3did, negatome_file).domainsplit_db

    if (params.negative_ppi_parquet != null) {
        domainsplit_db = INSERT_PPI_NEGATIVE_DDIS(
            domainsplit_db,
            file(params.negative_ppi_parquet),
            uniprot_id_mapping,
            params.negative_ppi_min_n_tested,
            params.negative_ppi_source_label,
            params.negative_sampling_strategy,
        ).domainsplit_db
    }

    if (params.smoke_test_n_ddis != null) {
        domainsplit_db = SMOKE_FILTER(domainsplit_db, params.smoke_test_n_ddis).domainsplit_db
    }

    emit:
    domainsplit_db
}
