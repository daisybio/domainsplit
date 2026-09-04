/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT_DDI_DATA -- download and parse every DDI source. 3did goes straight
    into the pre-initialised Domainsplit SQLite; every other source is parsed to
    a normalized TSV and inserted afterwards. The 3did SQLite stays internal to
    this subworkflow.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Two classes of source, and the split between them is load-bearing:

      * 3did (and, later, ppi-splitting's sampled negatives) are *protected*:
        they own their domain pairs outright, because they are the population
        ppi-splitting partitions into train/validation/test.
      * single_domain_ppi / PPIDM / negatome are *external*: they are parsed to
        `pfam_a, pfam_b, negative, source` TSVs here and inserted last, by
        INSERT_EXTERNAL_SOURCES, which drops any pair a protected source already
        claimed. That drop is what makes the external test set strictly unseen.

    The parse/insert split also exists so the external families are known
    *early* -- they feed the union Pfam fetch -- while their rows land *late*.
    INSERT_EXTERNAL_SOURCES therefore lives in `workflows/domainsplit.nf`, after
    ingest; this subworkflow only emits the parsed `external_ddis`.

    Add a new DDI source by:
      1. Adding its download/parse module, emitting the normalized TSV
         (bin/external_ddi_tsv.py).
      2. Mixing its TSV into `external_ddis` below.
      3. Tagging its rows with a unique source string; merge/drop behaviour then
         follows automatically from bin/ddi_db_utils.py.
----------------------------------------------------------------------------*/

include { DOWNLOAD_3DID_SQLITE      } from '../../../modules/local/3did/main.nf'
include { DOWNLOAD_NEGATOME         } from '../../../modules/local/negatome/main.nf'
include { INSERT_3DID               } from '../../../modules/local/insert_3did/main.nf'
include { BUILD_SWISSPROT_PFAM_MAP  } from '../../../modules/local/swissprot_map/main.nf'
include { PARSE_SINGLE_DOMAIN_PPI   } from '../../../modules/local/parse_single_domain_ppi/main.nf'
include { PARSE_PPIDM               } from '../../../modules/local/parse_ppidm/main.nf'
include { PARSE_NEGATOME            } from '../../../modules/local/parse_negatome/main.nf'
include { BUILD_CANDIDATE_NETWORK   } from '../../../modules/local/build_candidate_network/main.nf'

workflow COLLECT_DDI_DATA {
    take:
    domainsplit_db_in
    url_3did
    url_negatome
    swissprot_pfam_tsv
    hippie_tsv
    ppidm_tsv
    negative_ppi_parquet

    main:
    ch_versions = Channel.empty()

    if( !hippie_tsv || !ppidm_tsv || !negative_ppi_parquet ) {
        log.error "Required inputs missing: hippie_tsv, ppidm_tsv, and negative_ppi_parquet must be provided"
        exit 1
    }
    file_3did     = file(url_3did)
    sqlite_3did   = DOWNLOAD_3DID_SQLITE(file_3did).sqlite
    negatome_file = DOWNLOAD_NEGATOME(url_negatome).negatome

    // 1. 3did positives -- the only source written to the DB at this stage.
    domainsplit_db = INSERT_3DID(domainsplit_db_in, sqlite_3did).domainsplit_db

    // 2. external sources -> normalized TSVs (no DB writes)
    swissprot_map = BUILD_SWISSPROT_PFAM_MAP(swissprot_pfam_tsv).map
    single_domain_ddis = PARSE_SINGLE_DOMAIN_PPI(
        file(hippie_tsv),
        swissprot_map,
        params.hippie_min_score,
    ).ddis
    ppidm_ddis    = PARSE_PPIDM(file(ppidm_tsv), params.ppidm_classes).ddis
    negatome_ddis = PARSE_NEGATOME(negatome_file).ddis

    // Fixed order so the merged source lists and the logs are reproducible;
    // the merge itself is order-independent by construction.
    external_ddis = single_domain_ddis
        .combine(ppidm_ddis)
        .combine(negatome_ddis)
        .map { single, ppidm, negatome -> [single, ppidm, negatome] }

    // 3. candidate network of high-confidence non-interacting Pfam pairs, for
    //    ppi-splitting's `ilp_candidates` negative sampler. Built against the
    //    3did-only DB, so its family universe is exactly 3did's.
    candidates = BUILD_CANDIDATE_NETWORK(
        domainsplit_db,
        file(negative_ppi_parquet),
        params.negative_ppi_min_n_tested,
        params.negative_ppi_gene_mapping
            ? file(params.negative_ppi_gene_mapping, checkIfExists: true)
            : [],
    )

    ch_versions = ch_versions.mix(
        DOWNLOAD_3DID_SQLITE.out.versions,
        DOWNLOAD_NEGATOME.out.versions,
        INSERT_3DID.out.versions,
        BUILD_SWISSPROT_PFAM_MAP.out.versions,
        PARSE_SINGLE_DOMAIN_PPI.out.versions,
        PARSE_PPIDM.out.versions,
        PARSE_NEGATOME.out.versions,
        BUILD_CANDIDATE_NETWORK.out.versions,
    )

    emit:
    domainsplit_db
    external_ddis
    // What 3did offered before dedup. Emitted rather than published because
    // REPORT_DDI_ATTRITION is the only consumer and reads it from the channel.
    offered_counts    = INSERT_3DID.out.counts
    candidate_network = candidates.candidate_network
    pfam_mapping      = candidates.pfam_mapping
    versions          = ch_versions
}
