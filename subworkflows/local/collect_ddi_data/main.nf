/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT_DDI_DATA -- download and parse every DDI source, write into the
    pre-initialised Domainsplit SQLite. Downstream code consumes only the
    database; the 3did SQLite stays internal to this subworkflow.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Sources are inserted in a fixed order so that, on a duplicate domain pair,
    the earlier source wins the label (INSERT OR IGNORE):

        3did -> single-domain PPI -> PPIDM -> negatome
          -> [optional self-interaction removal]
          -> high-confidence non-PPI negatives (over 3did domains only)

    Add a new DDI source by:
      1. Adding its download/parse module.
      2. Slotting an INSERT_<source> call into the chain below (each collects its
         own unique Pfam IDs and bulk-creates missing domain rows via the shared
         bin/ddi_db_utils.py helper).
      3. Tagging its rows with a unique source string in domain_domain_interaction.
----------------------------------------------------------------------------*/

include { DOWNLOAD_3DID_SQLITE      } from '../../../modules/local/3did/main.nf'
include { DOWNLOAD_NEGATOME         } from '../../../modules/local/negatome/main.nf'
include { INSERT_3DID               } from '../../../modules/local/insert_3did/main.nf'
include { BUILD_SWISSPROT_PFAM_MAP  } from '../../../modules/local/swissprot_map/main.nf'
include { INSERT_SINGLE_DOMAIN_PPI  } from '../../../modules/local/insert_single_domain_ppi/main.nf'
include { INSERT_PPIDM              } from '../../../modules/local/insert_ppidm/main.nf'
include { INSERT_NEGATOME           } from '../../../modules/local/insert_negatome/main.nf'
include { REMOVE_SELF_INTERACTIONS  } from '../../../modules/local/remove_self_interactions/main.nf'
include { BUILD_PPI_NEGATIVE_POOL                            } from '../../../modules/local/build_ppi_negative_pool/main.nf'
include { SELECT_PPI_NEGATIVE_DANS as SELECT_DELETION        } from '../../../modules/local/select_ppi_negative_dans/main.nf'
include { SELECT_PPI_NEGATIVE_DANS as SELECT_RANDOM_ADDITION } from '../../../modules/local/select_ppi_negative_dans/main.nf'
include { INSERT_PPI_NEGATIVE_SELECTION                      } from '../../../modules/local/insert_ppi_negative_selection/main.nf'
include { SMOKE_FILTER              } from '../../../modules/local/smoke_filter/main.nf'

workflow COLLECT_DDI_DATA {
    take:
    domainsplit_db_in
    url_3did
    url_negatome
    url_uniprot_swissprot_pfam
    hippie_tsv
    ppidm_tsv
    negative_ppi_parquet

    main:
    if( !hippie_tsv || !ppidm_tsv || !negative_ppi_parquet ) {
        log.error "Required inputs missing: hippie_tsv, ppidm_tsv, and negative_ppi_parquet must be provided"
        exit 1
    }
    file_3did     = file(url_3did)
    sqlite_3did   = DOWNLOAD_3DID_SQLITE(file_3did).sqlite
    negatome_file = DOWNLOAD_NEGATOME(url_negatome).negatome

    // 1. 3did positives
    domainsplit_db = INSERT_3DID(domainsplit_db_in, sqlite_3did).domainsplit_db

    // 2-3. single-domain PPI positives (HIPPIE), using a reviewed-human SwissProt map
    swissprot_map = BUILD_SWISSPROT_PFAM_MAP(url_uniprot_swissprot_pfam).map
    domainsplit_db = INSERT_SINGLE_DOMAIN_PPI(
        domainsplit_db,
        file(hippie_tsv),
        swissprot_map,
        params.hippie_min_score,
    ).domainsplit_db

    // 4. PPIDM predicted positives (class kept as source)
    domainsplit_db = INSERT_PPIDM(
        domainsplit_db,
        file(ppidm_tsv),
        params.ppidm_classes,
    ).domainsplit_db

    // 5. negatome negatives
    domainsplit_db = INSERT_NEGATOME(domainsplit_db, negatome_file).domainsplit_db

    // 6. optional removal of all self-interactions
    if (!params.self_interaction) {
        domainsplit_db = REMOVE_SELF_INTERACTIONS(domainsplit_db).domainsplit_db
    }

    // 7. high-confidence non-PPI negatives via uncapped DANS (Cappelletti et al.
    //    vbae036), in two flavours that coexist under distinct source labels:
    //      * "deletion"        -- DANS over the PPI candidate pool, with the
    //                             positives reduced to the candidate-domain
    //                             universe (labels 3did_deletion /
    //                             inferred_ppi_screen_negative_for_deletion).
    //      * "random_addition" -- plain DANS over the full positive set (labels
    //                             3did_random_addition /
    //                             inferred_ppi_screen_negative_for_random_addition).
    //    The expensive, deterministic UniProt fetch + candidate-pool build runs
    //    once; each method is a single deterministic selection (no pick-best).
    pool = BUILD_PPI_NEGATIVE_POOL(
        domainsplit_db,
        file(negative_ppi_parquet),
        params.negative_ppi_min_n_tested,
        params.self_interaction,
    )

    del  = SELECT_DELETION('deletion', params.negative_ppi_seed, pool.neg_pool)
    rand = SELECT_RANDOM_ADDITION('random_addition', params.negative_ppi_seed, pool.neg_pool)

    inserted = INSERT_PPI_NEGATIVE_SELECTION(
        pool.domainsplit_db,
        pool.neg_pool,
        del.pairs,
        rand.pairs,
        del.score,
        rand.score,
    )
    domainsplit_db = inserted.domainsplit_db
    pfam_mapping   = pool.pfam_mapping

    if (params.smoke_test_n_ddis != null) {
        domainsplit_db = SMOKE_FILTER(domainsplit_db, params.smoke_test_n_ddis).domainsplit_db
    }

    emit:
    domainsplit_db
    pfam_mapping
}
