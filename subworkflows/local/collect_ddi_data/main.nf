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

// The source tokens `--external_test_sources` accepts. `PPIDM_<Class>` covers
// the classes parse_ppidm.py emits.
def externalSourceTokens() {
    return ['single_domain_ppi', 'PPIDM', 'PPIDM_Gold', 'PPIDM_Silver', 'PPIDM_Bronze', 'negatome']
}

// Whether this run needs an external source at all: the external sources feed
// only the external test set (and the conflict drop between themselves), so a
// source `--external_test_sources` does not select is neither downloaded nor
// parsed. Empty param = every source. `source` is one of single_domain_ppi,
// PPIDM, negatome.
def externalSourceWanted(source) {
    def tokens = params.external_test_sources.toString().tokenize(',').collect { t -> t.trim() }.findAll { t -> t }
    if (!tokens) {
        return true
    }
    def unknown = tokens.findAll { t -> !(t in externalSourceTokens()) }
    if (unknown) {
        error("--external_test_sources: unknown source ${unknown.join(', ')} (known: ${externalSourceTokens().join(', ')})")
    }
    if (source == 'PPIDM') {
        return tokens.any { t -> t == 'PPIDM' || t.startsWith('PPIDM_') }
    }
    return source in tokens
}

workflow COLLECT_DDI_DATA {
    take:
    domainsplit_db_in
    url_3did
    url_negatome
    swissprot_pfam_tsv
    hippie_tsv
    ppidm_tsv
    negative_ppi_parquet
    want_candidates   // Boolean: some split row asks for the `ilp_candidates` negative set

    main:
    ch_versions = Channel.empty()

    def want_single   = externalSourceWanted('single_domain_ppi')
    def want_ppidm    = externalSourceWanted('PPIDM')
    def want_negatome = externalSourceWanted('negatome')
    if (want_candidates && !negative_ppi_parquet) {
        error("--negative_ppi_parquet is required while a split row samples `ilp_candidates` negatives (see --ppi_splitting_multi_negset)")
    }
    if (want_single && !hippie_tsv) {
        error("--hippie_tsv is required while the external test set uses single_domain_ppi (see --external_test_sources)")
    }
    if (want_ppidm && !ppidm_tsv) {
        error("--ppidm_tsv is required while the external test set uses PPIDM (see --external_test_sources)")
    }
    file_3did     = file(url_3did)
    sqlite_3did   = DOWNLOAD_3DID_SQLITE(file_3did).sqlite

    // 1. 3did positives -- the only source written to the DB at this stage.
    domainsplit_db = INSERT_3DID(domainsplit_db_in, sqlite_3did).domainsplit_db

    // 2. external sources -> normalized TSVs (no DB writes)
    //    Only the sources --external_test_sources selects are fetched/parsed.
    def parsed = []
    if (want_single) {
        swissprot_map = BUILD_SWISSPROT_PFAM_MAP(swissprot_pfam_tsv).map
        parsed << PARSE_SINGLE_DOMAIN_PPI(
            file(hippie_tsv),
            swissprot_map,
            params.hippie_min_score,
        ).ddis
        ch_versions = ch_versions.mix(BUILD_SWISSPROT_PFAM_MAP.out.versions, PARSE_SINGLE_DOMAIN_PPI.out.versions)
    }
    if (want_ppidm) {
        parsed << PARSE_PPIDM(file(ppidm_tsv), params.ppidm_classes).ddis
        ch_versions = ch_versions.mix(PARSE_PPIDM.out.versions)
    }
    if (want_negatome) {
        parsed << PARSE_NEGATOME(DOWNLOAD_NEGATOME(url_negatome).negatome).ddis
        ch_versions = ch_versions.mix(DOWNLOAD_NEGATOME.out.versions, PARSE_NEGATOME.out.versions)
    }

    // Fixed order (single_domain_ppi, PPIDM, negatome) so the merged source
    // lists and the logs are reproducible; the merge itself is
    // order-independent by construction. combine() flattens, so the result is
    // one list of 1-3 TSVs.
    external_ddis = parsed.drop(1)
        .inject(parsed[0]) { acc, ch -> acc.combine(ch) }
        .map { tsvs -> (tsvs instanceof List) ? tsvs : [tsvs] }

    // 3. candidate network of high-confidence non-interacting Pfam pairs, for
    //    ppi-splitting's `ilp_candidates` negative sampler. Built against the
    //    3did-only DB, so its family universe is exactly 3did's. Only when some
    //    row samples from it: unmapped genes go to the UniProt REST API.
    if (want_candidates) {
        candidate_network = BUILD_CANDIDATE_NETWORK(
            domainsplit_db,
            file(negative_ppi_parquet),
            params.negative_ppi_min_n_tested,
            params.negative_ppi_gene_mapping
                ? file(params.negative_ppi_gene_mapping, checkIfExists: true)
                : [],
        ).candidate_network
        ch_versions = ch_versions.mix(BUILD_CANDIDATE_NETWORK.out.versions)
    }
    else {
        // A placeholder value, never staged: the caller hands the network only
        // to rows that want `ilp_candidates`, and here none does.
        candidate_network = channel.value('none')
    }

    ch_versions = ch_versions.mix(
        DOWNLOAD_3DID_SQLITE.out.versions,
        INSERT_3DID.out.versions,
    )

    emit:
    domainsplit_db
    external_ddis
    // What 3did offered before dedup. Emitted rather than published because
    // REPORT_DDI_ATTRITION is the only consumer and reads it from the channel.
    offered_counts    = INSERT_3DID.out.counts
    candidate_network
    versions          = ch_versions
}
