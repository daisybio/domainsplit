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
      * single_domain_ppi / PPIDM / INVERSE_PPIDM / negatome are *external*: they are parsed to
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
include { PARSE_SINGLE_DOMAIN_PPI   } from '../../../modules/local/parse_single_domain_ppi/main.nf'
include { PARSE_PPIDM               } from '../../../modules/local/parse_ppidm/main.nf'
include { PARSE_PPIDM as PARSE_INVERSE_PPIDM } from '../../../modules/local/parse_ppidm/main.nf'
include { EXTRACT_PPI_MITAB; BUILD_PPIDM_SOURCES; COMPUTE_PPIDM } from '../../../modules/local/ppidm/main.nf'
include { PARSE_NEGATOME            } from '../../../modules/local/parse_negatome/main.nf'
include { BUILD_CANDIDATE_NETWORK   } from '../../../modules/local/build_candidate_network/main.nf'

// The source tokens `--external_test_sources` accepts. `PPIDM_<Class>` covers
// the classes parse_ppidm.py emits; `INVERSE_PPIDM*` are the computed negatives
// (no Bronze: the lower tail has no "above threshold, not significant" class).
def externalSourceTokens() {
    return ['single_domain_ppi', 'PPIDM', 'PPIDM_Gold', 'PPIDM_Silver', 'PPIDM_Bronze',
            'INVERSE_PPIDM', 'INVERSE_PPIDM_Gold', 'INVERSE_PPIDM_Silver', 'negatome']
}

def externalTestTokens() {
    return params.external_test_sources.toString().tokenize(',').collect { t -> t.trim() }.findAll { t -> t }
}

// Whether this run needs an external source at all: the external sources feed
// only the external test set (and the conflict drop between themselves), so a
// source `--external_test_sources` does not select is neither downloaded nor
// parsed. Empty param = every source. `source` is one of single_domain_ppi,
// PPIDM, INVERSE_PPIDM, negatome.
//
// A legacy `--ppidm_tsv` replaces the computed positives, and computed
// negatives next to externally computed positives would mix two methods'
// label sets -- so INVERSE_PPIDM is refused when named explicitly and skipped
// when only implied by an empty token list.
def externalSourceWanted(source) {
    def tokens = externalTestTokens()
    def unknown = tokens.findAll { t -> !(t in externalSourceTokens()) }
    if (unknown) {
        error("--external_test_sources: unknown source ${unknown.join(', ')} (known: ${externalSourceTokens().join(', ')})")
    }
    if (source == 'INVERSE_PPIDM') {
        def named = tokens.any { t -> t == 'INVERSE_PPIDM' || t.startsWith('INVERSE_PPIDM_') }
        if (named && params.ppidm_tsv) {
            error("--external_test_sources names INVERSE_PPIDM, but --ppidm_tsv is set. The legacy --ppidm_tsv replaces the computed PPIDM positives, and its labels come from a different computation than the in-pipeline INVERSE_PPIDM negatives. Drop --ppidm_tsv to compute both, or drop the INVERSE_PPIDM tokens.")
        }
        return tokens ? named : !params.ppidm_tsv
    }
    if (!tokens) {
        return true
    }
    if (source == 'PPIDM') {
        return tokens.any { t -> t == 'PPIDM' || t.startsWith('PPIDM_') }
    }
    return source in tokens
}

// The classes PARSE_PPIDM (prefix `PPIDM`) or PARSE_INVERSE_PPIDM (prefix
// `INVERSE_PPIDM`) emits. When every selected token is a PPIDM/INVERSE_PPIDM
// one, explicit `<prefix>_<Class>` tokens narrow --ppidm_classes to the classes
// named; a class nothing selects would only widen the master DB and the union
// Pfam fetch / embedding workload. The two prefixes are disjoint by
// construction (compute_ppidm.py asserts it), so narrowing loses no conflict
// between them. Next to negatome or single_domain_ppi every configured class is
// parsed: INSERT_EXTERNAL_SOURCES drops a pair whose external labels disagree,
// and an unparsed class would hide that disagreement. A bare `<prefix>` token
// or an empty selection keeps --ppidm_classes as configured. Classes compare
// case-insensitively, as parse_ppidm.py does.
def ppidmParseClasses(prefix) {
    def tokens = externalTestTokens()
    def configured = params.ppidm_classes.toString().tokenize(',').collect { c -> c.trim() }.findAll { c -> c }
    def picked = tokens.findAll { t -> t.startsWith(prefix + '_') }.collect { t -> t.substring(prefix.length() + 1) }
    def missing = picked.findAll { c -> !configured.any { k -> k.equalsIgnoreCase(c) } }
    if (missing) {
        error("--external_test_sources selects ${missing.collect { c -> prefix + '_' + c }.join(', ')}, but --ppidm_classes (${params.ppidm_classes}) does not parse that class")
    }
    def ppidm_only = tokens.every { t -> t.startsWith('PPIDM') || t.startsWith('INVERSE_PPIDM') }
    def classes = !tokens || prefix in tokens || !ppidm_only
        ? params.ppidm_classes.toString()
        : configured.findAll { k -> picked.any { c -> c.equalsIgnoreCase(k) } }.join(',')
    log.info "PARSE_${prefix} classes: ${classes}"
    return classes
}

// Under `--external_test_balance uniform` the larger label is subsampled to the
// smaller, so an explicit selection naming only one label would fail late in
// BUILD_EXTERNAL_TEST after the whole build. Refuse it at startup instead
// (build_external_test.py keeps its own check for a label lost to filtering).
def checkExternalTestBalance() {
    def tokens = externalTestTokens()
    if (!tokens || params.external_test_balance.toString() != 'uniform') {
        return
    }
    def negative = tokens.findAll { t -> t == 'negatome' || t.startsWith('INVERSE_PPIDM') }
    if (!negative || negative.size() == tokens.size()) {
        error("--external_test_balance uniform needs both labels, but --external_test_sources (${tokens.join(',')}) selects only ${negative ? 'negatives' : 'positives'}. Add a ${negative ? 'positive (single_domain_ppi, PPIDM*)' : 'negative (INVERSE_PPIDM*, negatome)'} source, or pass --external_test_balance none.")
    }
}

// Whether COMPUTE_PPIDM runs: computed positives are wanted (no legacy
// --ppidm_tsv) or computed negatives are.
def ppidmComputed() {
    return (externalSourceWanted('PPIDM') && !params.ppidm_tsv) || externalSourceWanted('INVERSE_PPIDM')
}

// STRING's per-organism detailed links files for `ppidm_string_taxa`, from
// `ppidm_string_url_template` (`{taxid}`, `{release}` placeholders; a local
// path template works the same way).
def ppidmStringFiles() {
    def taxa = params.ppidm_string_taxa.toString().tokenize(',').collect { t -> t.trim() }.findAll { t -> t }
    if (!taxa) {
        error("--ppidm_string_taxa is empty, so the two STRING sources of PPIDM would have no PPIs. Give at least one NCBI taxon id (default: 9606,4932,10090,7227,6239,3702,511145).")
    }
    return taxa.collect { t ->
        params.ppidm_string_url_template.toString()
            .replace('{taxid}', t)
            .replace('{release}', params.string_release.toString())
    }
}

// Every remote/manual file the PPIDM computation reads, for the manifest.
def ppidmInputUrls() {
    def dats = params.ppidm_uniprot_dat.toString().tokenize(',').collect { u -> u.trim() }.findAll { u -> u }
    return dats + [params.url_biogrid_mitab, params.url_intact_mitab, params.url_hippie,
                   params.url_sifts_chain_uniprot, params.url_3did_flat] + ppidmStringFiles()
}

workflow COLLECT_DDI_DATA {
    take:
    domainsplit_db_in
    url_3did
    url_negatome
    swissprot_pfam_tsv
    hippie        // --url_hippie (URL or path); read by single_domain_ppi and PPIDM
    ppidm_tsv     // legacy precomputed PPIDM positives, or null to compute them
    negative_ppi_parquet
    want_candidates   // Boolean: some split row asks for the `ilp_candidates` negative set

    main:
    ch_versions = Channel.empty()

    def want_single   = externalSourceWanted('single_domain_ppi')
    def want_ppidm    = externalSourceWanted('PPIDM')
    def want_inverse  = externalSourceWanted('INVERSE_PPIDM')
    def want_negatome = externalSourceWanted('negatome')
    def compute_ppidm = ppidmComputed()
    checkExternalTestBalance()
    if (want_candidates && !negative_ppi_parquet) {
        error("--negative_ppi_parquet is required while a split row samples `ilp_candidates` negatives (see --ppi_splitting_multi_negset)")
    }
    if ((want_single || compute_ppidm) && !hippie) {
        error("--url_hippie is required while the external test set uses single_domain_ppi or a computed PPIDM/INVERSE_PPIDM source (see --external_test_sources)")
    }
    file_3did     = file(url_3did)
    sqlite_3did   = DOWNLOAD_3DID_SQLITE(file_3did).sqlite

    // 1. 3did positives -- the only source written to the DB at this stage.
    domainsplit_db = INSERT_3DID(domainsplit_db_in, sqlite_3did).domainsplit_db

    // 1b. PPIDM / INVERSE_PPIDM, computed from raw PPI downloads against the
    //     3did-only database as gold standard. Labels are independent of
    //     `instance_tier`: the accession->Pfam map is always all-species
    //     reviewed UniProt (`ppidm_uniprot_dat`), not PARSE_SWISSPROT's
    //     tier-filtered one.
    ppidm_report  = channel.empty()
    ppidm_sources = channel.empty()
    if (compute_ppidm) {
        mitab = EXTRACT_PPI_MITAB(channel.of(
            ['biogrid', file(params.url_biogrid_mitab)],
            ['intact',  file(params.url_intact_mitab)],
        )).ppis
        BUILD_PPIDM_SOURCES(
            params.ppidm_uniprot_dat.toString().tokenize(',').collect { u -> file(u.trim()) },
            mitab.filter { fmt, f -> fmt == 'biogrid' }.map { fmt, f -> f },
            mitab.filter { fmt, f -> fmt == 'intact' }.map { fmt, f -> f },
            file(hippie),
            ppidmStringFiles().collect { u -> file(u) },
            file(params.url_sifts_chain_uniprot),
            file(params.url_3did_flat),
        )
        COMPUTE_PPIDM(
            BUILD_PPIDM_SOURCES.out.sources,
            BUILD_PPIDM_SOURCES.out.protein_pfam,
            domainsplit_db,
        )
        ppidm_report  = COMPUTE_PPIDM.out.report
        ppidm_sources = BUILD_PPIDM_SOURCES.out.report
        ch_versions = ch_versions.mix(
            EXTRACT_PPI_MITAB.out.versions.first(),
            BUILD_PPIDM_SOURCES.out.versions,
            COMPUTE_PPIDM.out.versions,
        )
    }

    // 2. external sources -> normalized TSVs (no DB writes)
    //    Only the sources --external_test_sources selects are fetched/parsed.
    def parsed = []
    if (want_single) {
        parsed << PARSE_SINGLE_DOMAIN_PPI(
            file(hippie),
            swissprot_pfam_tsv,
            params.hippie_min_score,
        ).ddis
        ch_versions = ch_versions.mix(PARSE_SINGLE_DOMAIN_PPI.out.versions)
    }
    if (want_ppidm) {
        parsed << PARSE_PPIDM(
            ppidm_tsv ? file(ppidm_tsv, checkIfExists: true) : COMPUTE_PPIDM.out.ppidm,
            ppidmParseClasses('PPIDM'),
        ).ddis
        ch_versions = ch_versions.mix(PARSE_PPIDM.out.versions)
    }
    if (want_inverse) {
        parsed << PARSE_INVERSE_PPIDM(COMPUTE_PPIDM.out.inverse_ppidm, ppidmParseClasses('INVERSE_PPIDM')).ddis
        ch_versions = ch_versions.mix(PARSE_INVERSE_PPIDM.out.versions)
    }
    if (want_negatome) {
        // A local path is staged (so -resume sees its content); a URL keeps the
        // download process, whose cert fallback the default host needs (missing intermediate).
        def remote_negatome = url_negatome.toString() ==~ /^(https?|ftp):\/\/.*/
        parsed << PARSE_NEGATOME(
            remote_negatome ? DOWNLOAD_NEGATOME(url_negatome).negatome : file(url_negatome, checkIfExists: true)
        ).ddis
        ch_versions = ch_versions.mix(PARSE_NEGATOME.out.versions)
        if (remote_negatome) {
            ch_versions = ch_versions.mix(DOWNLOAD_NEGATOME.out.versions)
        }
    }

    // Fixed order (single_domain_ppi, PPIDM, INVERSE_PPIDM, negatome) so the merged source
    // lists and the logs are reproducible; the merge itself is
    // order-independent by construction. combine() flattens, so the result is
    // one list of 1-4 TSVs.
    external_ddis = parsed.drop(1)
        .inject(parsed[0]) { acc, ch -> acc.combine(ch) }
        .map { tsvs -> (tsvs instanceof List) ? tsvs : [tsvs] }

    // 3. candidate network of high-replication non-interacting Pfam pairs, for
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
        candidate_pfam_mapping = BUILD_CANDIDATE_NETWORK.out.pfam_mapping
        ch_versions = ch_versions.mix(BUILD_CANDIDATE_NETWORK.out.versions)
    }
    else {
        // A placeholder value, never staged: the caller hands the network only
        // to rows that want `ilp_candidates`, and here none does.
        candidate_network = channel.value('none')
        candidate_pfam_mapping = channel.value([])
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
    // The gene->Pfam JSON the network was built from (live UniProt unless pinned),
    // for the manifest's SHA-256; [] when no row samples `ilp_candidates`.
    candidate_pfam_mapping
    // COMPUTE_PPIDM's report and BUILD_PPIDM_SOURCES' per-source counts, for
    // the manifest; empty channels when PPIDM is not computed.
    ppidm_report
    ppidm_sources
    versions          = ch_versions
}
