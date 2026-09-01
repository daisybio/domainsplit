/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_domainsplit_pipeline'
include { INIT_DOMAINSPLIT_DB         } from '../modules/local/init_domainsplit_db/main.nf'
include { PARSE_SWISSPROT             } from '../modules/local/parse_swissprot/main.nf'
include { COLLECT_DDI_DATA            } from '../subworkflows/local/collect_ddi_data/main.nf'
include { EXPORT_UNION_FAMILIES       } from '../modules/local/export_union_families/main.nf'
include { EXPORT_SPLIT_DDIS           } from '../modules/local/export_split_ddis/main.nf'
include { INGEST_SPLITS               } from '../subworkflows/local/ingest_splits/main.nf'
include { INSERT_EXTERNAL_SOURCES     } from '../modules/local/insert_external_sources/main.nf'
include { PRUNE_UNREPRESENTED_DDIS    } from '../modules/local/prune_unrepresented_ddis/main.nf'
include { FETCH_STRING_LINKS          } from '../modules/local/fetch_string_links/main.nf'
include { REPORT_DDI_ATTRITION        } from '../modules/local/report_ddi_attrition/main.nf'
include { generate_domain_embeddings  } from '../modules/local/domain_embeddings/main.nf'
include { ENRICH_DDI_DATABASE         } from '../subworkflows/local/enrich_ddi_database/main.nf'
include { BUILD_EXTERNAL_TEST         } from '../modules/local/build_external_test/main.nf'
include { SUBSET_SPLIT_DB             } from '../modules/local/subset_split_db/main.nf'

// The imported splitting pipeline (git submodule, pinned in .gitmodules), plus
// the three data-prep processes we drive ourselves so the Pfam fetch happens
// exactly once, over the union of every family this run touches.
include { PPI_SPLITTING               } from '../subworkflows/external/ppi-splitting/main.nf'
include { FETCH_DOMAIN_META           } from '../subworkflows/external/ppi-splitting/processes/data_prep.nf'
include { GET_LENGTHS                 } from '../subworkflows/external/ppi-splitting/processes/data_prep.nf'
include { SUBSET_DOMAIN_DATA          } from '../subworkflows/external/ppi-splitting/processes/data_prep.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SPLITTING STRATEGIES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    One row per positive split. `negsets` maps ppi-splitting's negative-set name
    to the output directory the resulting databases land in -- that mapping lives
    here and nowhere else, because only this file knows which dataset row asked
    for which negative set. Every row splits the *same* 3did positives; the rows
    differ in how, and in what the negatives look like.

      random           naive baseline: random split, uniform negatives.
      minimal_leakage  ILP split, ILP negatives; the main product.
      external_test    same as minimal_leakage but with no internal test set --
                       its test split is BUILD_EXTERNAL_TEST's external DDIs,
                       which are unseen by construction (see §2 of the plan).

    `ilp_candidates` is the ILP sampler restricted to the candidate network of
    experimentally-tested-but-not-interacting Pfam pairs, so it yields
    high-confidence negative instances -- hence the `_hcni` directories. It needs
    ppi-splitting's comma-separated `negative_sampling_method` fan-out, which the
    pinned submodule has (`parseNegsets` in its main.nf), so both ILP rows produce
    two negative sets on one positive split and the run yields 5 method
    directories / 18 split databases.
----------------------------------------------------------------------------*/

// A function, not a script-level `def` variable: those are locals of the script's
// own run method and are not reliably visible from inside a workflow body.
def splitRows() {
    return [
        [
            id            : 'random',
            split_method  : 'random',
            negsets       : ['uniform': 'random'],
            internal_test : true,
        ],
        [
            id            : 'minimal_leakage',
            split_method  : 'ilp',
            negsets       : ['ilp': 'minimal_leakage', 'ilp_candidates': 'minimal_leakage_hcni'],
            internal_test : true,
        ],
        [
            id            : 'external_test',
            split_method  : 'ilp',
            negsets       : ['ilp': 'external_test', 'ilp_candidates': 'external_test_hcni'],
            internal_test : false,
        ],
    ]
}

// The negative sets a row actually produces on this run. Normally all of them:
// `params.ppi_splitting_multi_negset` exists so a `--split_only` debug run can
// drop back to one, because ppi-splitting rejects `--split_only` unless
// `negative_sampling_method` is exactly "ilp" (its main.nf). With the param off,
// each row keeps only its *first* negative set and the `_hcni` directories are
// not produced.
def activeNegsets(row) {
    // `.toString().toBoolean()`, not the raw param: a CLI `--ppi_splitting_multi_negset
    // false` arrives as the *String* "false", and every non-empty String is truthy in
    // Groovy -- so testing the param directly would silently keep both negative sets on
    // exactly the run that asked for one. It mattered less when the default was `false`
    // (turning it *on* passes "true", which is truthy by luck); now that the default is
    // `true`, turning it off is the operation that has to work.
    def multi = params.ppi_splitting_multi_negset.toString().toBoolean()
    return multi ? row.negsets.keySet() as List : [row.negsets.keySet().first()]
}

// ppi-splitting labels its splits train/val/test_balanced/test_realistic and
// the ingest scripts rename val -> validation. A row with no internal test set
// contributes train/validation here and gets its `test` from BUILD_EXTERNAL_TEST.
def splitsOf(row) {
    return row.internal_test ? ['train', 'validation', 'test_balanced', 'test_realistic'] : ['train', 'validation']
}

// Every per-dataset override ppi-splitting reads out of `meta`. Mirrors
// buildDatasetsChannel() in its main.nf: the keys are fixed, and mutating the
// map downstream would break every join() in there, all of which key on the
// whole map.
def metaOf(row) {
    def negsets = activeNegsets(row)
    return [
        id                       : row.id,
        embedding_model          : params.embedding_model,
        cdhit_identity           : params.cdhit_identity,
        cdhit_wordsize           : params.cdhit_wordsize,
        split_method             : row.split_method,
        edge_weight              : params.edge_weight,
        kahip_k                  : params.kahip_k,
        ilp_kahip_k              : params.ilp_kahip_k,
        train_split              : row.internal_test ? params.split_train_fraction : (1.0 - params.split_val_fraction),
        val_split                : params.split_val_fraction,
        test_split               : row.internal_test ? params.split_test_fraction : 0.0,
        ilp_epsilon              : params.ilp_epsilon,
        ilp_max_sec              : params.ilp_max_sec,
        negative_sampling_method : negsets.join(','),
        neg_ilp_time_limit       : params.neg_ilp_time_limit,
        neg_ilp_lambda_degree    : params.neg_ilp_lambda_degree,
        neg_ilp_lambda_taxon_pair: params.neg_ilp_lambda_taxon_pair,
        neg_ilp_lambda_self_loop : params.neg_ilp_lambda_self_loop,
        // GO terms describe proteins, not domain families, so the Jaccard bias
        // term has nothing to match on in DDI mode; ppi-splitting forces it to 0
        // itself, and restating it here keeps the two metas identical.
        neg_ilp_lambda_jaccard   : 0,
    ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow DOMAINSPLIT {
main:
    ch_versions = Channel.empty()

    //
    // One knob, three derived values, and every failure mode here is silence.
    //
    // `instance_tier` is the protein universe. It decides which UniProt flat
    // files are downloaded (`uniprot_dat_urls`), which taxa survive parsing
    // (`swissprot_taxon_ids`) and which strata FETCH_DOMAIN_META may sample
    // (`instance_tiers`, ppi-splitting's vocabulary). Those three are derived in
    // nextflow.config -- not here -- because FETCH_DOMAIN_META lives in the
    // read-only submodule and reads `params.instance_tiers` straight out of its
    // own script block, and `params` is read-only by the time a workflow body
    // runs.
    //
    // What is left here is the checking. Each guard exists because the failure it
    // catches otherwise produces a green run with fewer families and nothing
    // saying why. Checked before the first task is submitted, so it costs nothing.
    //
    def tier = params.instance_tier.toString()

    if (tier == 'all_species_any_review_status') {
        error("""--instance_tier all_species_any_review_status is not implemented.

It needs the full TrEMBL flat file (uniprot_trembl.dat.gz, 110 GB, ~250M entries),
which FETCH_DOMAIN_META cannot hold: it parses the protein universe into an in-memory
dict. Implementing it requires a cached accession -> (taxon, reviewed, byte offset)
index and offset-seek sequence retrieval.

See https://github.com/daisybio/domainsplit/issues/4
Use all_species_reviewed or human_any_review_status.""")
    }
    // The retired vocabulary, mapped rather than aliased. `any` used to mean "four
    // strata over a Swiss-Prot universe", i.e. the new `all_species_reviewed` --
    // NOT `all_species_any_review_status`. Aliasing it would silently change the
    // universe of a working invocation, so it is an error with the mapping instead.
    if (tier in ['human_only', 'any']) {
        error("--instance_tier '${tier}' is the retired vocabulary. 'human_only' is now 'human_reviewed'; 'any' is now 'all_species_reviewed' (it meant a Swiss-Prot universe over every species -- NOT 'all_species_any_review_status', which adds TrEMBL and is unimplemented).")
    }
    if (!params.instance_tiers) {
        error("--instance_tier '${tier}' is not a known protein universe. One of: human_reviewed, all_species_reviewed, human_any_review_status, all_species_any_review_status (unimplemented, see daisybio/domainsplit#4).")
    }

    // `--swissprot_taxon_ids ''` on the command line is parsed by Nextflow as a
    // bare flag and arrives as Boolean true, which would then stringify to "true"
    // and filter every entry out. There is no reason to pass it empty any more --
    // an `all_species_*` tier derives it -- so say that rather than let it through.
    if (params.swissprot_taxon_ids instanceof Boolean) {
        error("--swissprot_taxon_ids was given no value; Nextflow parsed it as a bare boolean flag. An empty taxon filter is what `--instance_tier all_species_reviewed` already derives, so drop the flag; to set it explicitly to empty, use a -params-file.")
    }

    def taxa = params.swissprot_taxon_ids.toString().trim()
    def uniprot_dats = (params.uniprot_dat_urls instanceof List ? params.uniprot_dat_urls : params.uniprot_dat_urls.toString().tokenize(','))
        .collect { u -> u.toString().trim() }
        .findAll { u -> u }

    if (!uniprot_dats) {
        error("--uniprot_dat_urls resolved to nothing for --instance_tier ${tier}. Pass a comma-separated list of UniProt flat files, Swiss-Prot first.")
    }
    if (tier.startsWith('human_') && taxa && !(taxa.tokenize(',').collect { t -> t.trim() }.contains('9606'))) {
        error("--instance_tier ${tier} samples human instances, but --swissprot_taxon_ids is '${taxa}', which does not include 9606 -- the universe would contain no human protein and every family would end with zero instances. Add 9606, or use an all_species_* tier.")
    }
    // An unreviewed stratum that can never fill looks exactly like a successful run
    // with fewer families. The authoritative check is upstream's, against the
    // *parsed* universe; this one exists only to catch the derived-default case
    // before a 4.7 GB regions stream, so it fires only when the file list is
    // *provably* reviewed-only -- every entry named `uniprot_sprot*`, UniProt's
    // own Swiss-Prot filenames. A locally-named file it cannot classify passes
    // here and is judged upstream on its contents, which is the check that can
    // actually be right.
    def reviewed_only = uniprot_dats.every { u -> u.toLowerCase() ==~ /.*uniprot_sprot.*/ }
    if (tier.endsWith('_any_review_status') && reviewed_only) {
        error("--instance_tier ${tier} asks for an unreviewed stratum, but every resolved --uniprot_dat_urls entry is a Swiss-Prot flat file:\n  ${uniprot_dats.join('\n  ')}\nAn unreviewed stratum that can never fill produces a green run with fewer families and nothing saying why. Add the matching uniprot_trembl_*.dat.gz.")
    }
    // Not an error: this is how a human+mouse run is expressed.
    if (tier.startsWith('all_species_') && taxa) {
        log.info "--instance_tier ${tier} with --swissprot_taxon_ids '${taxa}': the universe is deliberately narrowed to those taxa, not every species."
    }

    // The whole derivation on one line of the run log. Four values that must
    // agree, three of them invisible because they are derived -- so a run whose
    // family count surprises someone six months later can be diagnosed from its
    // log rather than from a guess about which defaults were in effect.
    log.info "protein universe: instance_tier=${tier} -> tiers=${params.instance_tiers}, taxa=${taxa ?: 'all species'}, ${uniprot_dats.size()} UniProt flat file(s): ${uniprot_dats.collect { u -> file(u).name }.join(', ')}"

    input_pfam2go = file(params.url_pfam2go)

    //
    // One pass over the UniProt flat files, six consumers. See the
    // `uniprot_dat_urls` comment in nextflow.config for why this is a parse of
    // static FTP files rather than two REST stream queries.
    //
    // The sixth consumer is FETCH_DOMAIN_META, which needs them as the protein
    // universe -- handed the raw files rather than a PARSE_SWISSPROT output,
    // because it needs sequences and taxa for every accession Pfam might name,
    // which is a superset of what the enrichment outputs carry.
    //
    PARSE_SWISSPROT(
        uniprot_dats.collect { u -> file(u) },
        params.swissprot_taxon_ids,
    )
    ch_versions = ch_versions.mix(PARSE_SWISSPROT.out.versions)

    input_uniprot_go_terms  = PARSE_SWISSPROT.out.go_terms
    input_uniprot_sequences = PARSE_SWISSPROT.out.sequences
    input_swissprot_pfam    = PARSE_SWISSPROT.out.pfam_tsv
    // Entry -> STRING id for every species the taxon filter admits, replacing the
    // per-organism idmapping download INSERT_PPI used to take.
    input_string_map        = PARSE_SWISSPROT.out.string_map

    empty_db = INIT_DOMAINSPLIT_DB().domainsplit_db

    //
    // 3did into the database; every other source parsed to a normalized TSV.
    //
    COLLECT_DDI_DATA(
        empty_db,
        params.url_3did,
        params.url_negatome,
        input_swissprot_pfam,
        params.hippie_tsv,
        params.ppidm_tsv,
        params.negative_ppi_parquet,
    )

    domainsplit_db_ddi = COLLECT_DDI_DATA.out.domainsplit_db
    external_ddis      = COLLECT_DDI_DATA.out.external_ddis

    //
    // ONE Pfam fetch, over the union of 3did and external-source families.
    //
    // External families need instances too -- BUILD_EXTERNAL_TEST represents
    // their DDIs as instance pairs, and domain_protein_map has to cover them --
    // but they must not enter the split population, so the fetched data is
    // subset back to 3did families before it reaches PPI_SPLITTING. Otherwise
    // BLAST would run all-vs-all over instances of families that are never split.
    //
    families = EXPORT_UNION_FAMILIES(domainsplit_db_ddi, external_ddis).families

    // Release-wide Pfam reference files and the download cache are run-global
    // params on ppi-splitting's side; resolve them the same way DATA_PREP_DDI
    // does, including absolutising the cache (FETCH_DOMAIN_META *writes* to it,
    // so a relative path would land in the task's work dir and be deleted).
    clans_ch = params.pfam_clans
        ? channel.value(file(params.pfam_clans, checkIfExists: true))
        : channel.value([])
    regions_ch = params.pfam_regions
        ? channel.value(file(params.pfam_regions, checkIfExists: true))
        : channel.value(file(params.url_pfam_regions))
    // The protein universe. The same flat files PARSE_SWISSPROT reads, handed over
    // raw: FETCH_DOMAIN_META needs a sequence, a taxon and a review flag for every
    // accession Pfam might name, which is wider than any of the parsed outputs.
    // A list, Swiss-Prot first -- the upstream parser takes the first writer on a
    // duplicate accession, and the reviewed record is the one to keep.
    sprot_ch = channel.value(uniprot_dats.collect { u -> file(u) })
    cache_dir = params.cache_dir ? file(params.cache_dir).toAbsolutePath().toString() : ''
    if (cache_dir && !file(cache_dir).exists()) {
        log.warn "--cache_dir ${cache_dir} does not exist yet; it will be created."
    }

    fetched = FETCH_DOMAIN_META(
        families.map { fams -> tuple([id: 'union'], fams) },
        clans_ch,
        regions_ch,
        sprot_ch,
        channel.value(cache_dir),
    )

    // The families Pfam had nothing usable for, one row per family with a
    // `reason` column. Taken from the channel, never read back from the
    // published path. Three reasons and they are not interchangeable:
    // `no_eligible_instances` is a consequence of `instance_tier`, which defaults
    // to `human_reviewed`, so it can fire in bulk and would bury the
    // other two -- `dead` and `not_in_pfam` are facts about Pfam, and a
    // mistyped accession arrives as `not_in_pfam`. So they are counted apart
    // and only the Pfam-fact reasons warn.
    //
    // Note this is *our* fetch's report, not PPI_SPLITTING.out.dropped_families:
    // that one carries the synthetic [id: '_shared'] meta (which matches no
    // dataset, so join()ing it against a per-dataset channel silently drops the
    // branch) and is empty here anyway, because every dataset row supplies
    // precomputed `domain_instances` and DATA_PREP_DDI never fetches. Nothing
    // may block on either of them: a run that drops nothing still emits the
    // header-only file.
    dropped_families = fetched.dropped_families.map { _meta, f -> f }
    dropped_families
        .map { f -> f.readLines().drop(1).findAll { line -> line.trim() }.countBy { line -> line.tokenize('\t')[-1] } }
        .subscribe { counts ->
            def total = counts.values().sum() ?: 0
            if (total) {
                log.info "FETCH_DOMAIN_META dropped ${total} families: ${counts.sort().collect { r, n -> "${r}=${n}" }.join(', ')} (see dropped_families.tsv; their DDIs leave the database in PRUNE_UNREPRESENTED_DDIS)"
            }
            def pfam_facts = counts.findAll { reason, _n -> reason in ['dead', 'not_in_pfam'] }
            if (pfam_facts) {
                log.warn "${pfam_facts.collect { r, n -> "${n} families ${r}" }.join(', ')} -- these are not an instance_tier effect. Check for mistyped or retired accessions in the requested family set."
            }
        }

    union_sequences = fetched.sequences.map      { _meta, f -> f }
    union_species   = fetched.species.map        { _meta, f -> f }
    union_instances = fetched.instances.map      { _meta, f -> f }
    union_go        = fetched.go_annotations.map { _meta, f -> f }
    union_lengths   = GET_LENGTHS(fetched.sequences).map { _meta, f -> f }

    // The split population: 3did positives only, restricted to families that
    // actually resolved to instances. Under `instance_tier = human_reviewed` a
    // family with no human Swiss-Prot region ends with none, and its DDIs leave
    // the database later, in PRUNE_UNREPRESENTED_DDIS.
    split_ddis = EXPORT_SPLIT_DDIS(domainsplit_db_ddi, union_instances, '3did').ddis

    subset = SUBSET_DOMAIN_DATA(
        split_ddis.map { ddis -> tuple([id: '3did'], ddis) },
        union_sequences,
        union_go,
        union_species,
        union_lengths,
        union_instances,
    )

    //
    // Splitting + negative sampling, three dataset rows over one shared
    // positive set and one shared precomputed data prep.
    //
    candidate_network = COLLECT_DDI_DATA.out.candidate_network

    datasets_ch = channel.fromList(splitRows())
        .combine(split_ddis)
        .combine(subset.sequences.map { _meta, f -> f })
        .combine(subset.species.map   { _meta, f -> f })
        .combine(subset.instances.map { _meta, f -> f })
        .combine(candidate_network)
        .map { row, ppis, sequences, species, instances, candidates ->
            // A candidate network is passed only to rows that actually ask for
            // the `ilp_candidates` negative set; supplying it otherwise makes
            // ppi-splitting warn and ignore it. An absent optional file is [],
            // never null -- null blows up staging.
            def wants_candidates = 'ilp_candidates' in activeNegsets(row)
            def files = [
                ppis             : ppis,
                sequences        : sequences,
                go_annotations   : [],
                species          : species,
                domain_instances : instances,
                blast_results    : [],
                candidate_network: wants_candidates ? candidates : [],
                partition        : [],
                node_mapping     : [],
            ]
            tuple(metaOf(row), files)
        }

    ppi = PPI_SPLITTING(datasets_ch)

    // (meta, negset, label, csv) -> (method, split, csv). `test_balanced` /
    // `test_realistic` of the external_test row are empty files (its test
    // fraction is 0) and are dropped here.
    def negsetMethod = splitRows().collectEntries { row -> [(row.id): row.negsets] }
    def internalTest = splitRows().collectEntries { row -> [(row.id): row.internal_test] }

    family_splits = ppi.labelled
        .filter { meta, _negset, label, _csv -> internalTest[meta.id] || !label.startsWith('test') }
        .map    { meta, negset, label, csv -> tuple(negsetMethod[meta.id][negset], label, csv) }
    instance_splits = ppi.labelled_inst
        .filter { meta, _negset, label, _csv -> internalTest[meta.id] || !label.startsWith('test') }
        .map    { meta, negset, label, csv -> tuple(negsetMethod[meta.id][negset], label, csv) }

    //
    // Ingest, then the external sources, then prune -- in that order.
    //
    // The order is the insert-rule contract in bin/ddi_db_utils.py: 3did and the
    // sampled negatives own their pairs, so they must be in the database before
    // any external row is offered. INGEST_SAMPLED_NEGATIVES fails the run if it
    // finds an external source already holding a pair, so a mis-wire here is
    // loud rather than silent.
    //
    ingested = INGEST_SPLITS(
        domainsplit_db_ddi,
        union_instances,
        union_sequences,
        family_splits,
        instance_splits,
    )

    inserted = INSERT_EXTERNAL_SOURCES(ingested.domainsplit_db, external_ddis)
    pruned   = PRUNE_UNREPRESENTED_DDIS(inserted.domainsplit_db)

    domainsplit_db = pruned.domainsplit_db

    //
    // STRING links for every species in the database.
    //
    // Placed here because the `protein` table is final by this point -- ingest
    // created the rows, the prune removed the DDIs but not the proteins -- and
    // ENRICH_DDI_DATABASE, the consumer, runs after BUILD_EXTERNAL_TEST below.
    //
    // The taxon list is derived from the STRING ids themselves, not from
    // `protein.taxon_id`: a STRING id is `9606.ENSP00000269305` and its prefix is
    // STRING's own species-level taxid, while UniProt's OX can be a strain-level
    // id with no STRING file behind it. Taking the prefix makes the join exact by
    // construction.
    //
    // `--url_string` short-circuits the whole thing with one explicit links file,
    // which is what `-profile test` uses to stay offline.
    //
    if (params.url_string) {
        input_string = channel.value(file(params.url_string))
        log.info "--url_string is set: using ${params.url_string} verbatim; FETCH_STRING_LINKS is skipped and only the species that file covers are enriched."
    }
    else {
        input_string = FETCH_STRING_LINKS(
            domainsplit_db,
            input_string_map,
            channel.value(cache_dir),
        ).links
        ch_versions = ch_versions.mix(FETCH_STRING_LINKS.out.versions)
    }

    //
    // Where each source's DDIs went, in one table.
    //
    // PPI_SPLITTING's own attrition waterfall begins at `split_ddis` -- so the
    // 3did DDIs that `instance_tier` stranded before the splitting
    // stage are, correctly, absent from it, and there was no published artifact
    // showing them at all. The counters existed only in three `.command.log`
    // files, which do not survive the work directory.
    //
    // Keyed on the pruned DB, so `surviving` means "in the published master".
    // Every other input is a file that is published anyway, which is what makes
    // the report reproducible from a finished run.
    //
    attrition = REPORT_DDI_ATTRITION(
        domainsplit_db,
        pruned.pruned,
        split_ddis,
        external_ddis,
        inserted.conflicts,
        COLLECT_DDI_DATA.out.offered_counts,
        '3did',
    )

    //
    // The external test set: DDIs from PPIDM / single_domain_ppi / negatome,
    // none of which can have appeared in any split of any set. Written as the
    // `test` split of every external_test method directory.
    //
    // Deliberately *before* enrichment. BUILD_EXTERNAL_TEST clones the master,
    // and enrichment only ever grows it -- protein sequences, GO terms and the
    // STRING network -- so cloning first is the smaller copy. (It used to be a
    // far larger win: enrichment wrote per-residue embedding blobs that were
    // ~99% of the file. Those are gone, published as HDF5 instead, and the
    // ordering is now merely the cheaper one rather than the only tolerable one.)
    //
    // Safe because the two write disjoint tables and the read sets do not move:
    // this step reads `domain_domain_interaction`, `domain` and
    // `domain_protein_map` WHERE instance_id IS NOT NULL, and writes only
    // `ddi_split_membership`. Enrichment writes `protein`, `protein_go_terms`,
    // `domain_go_terms` and `protein_protein_interaction` and now touches
    // `domain_protein_map` not at all -- INSERT_DOMAIN_PROTEIN_MAPPING went with
    // the embedding columns. Sampling is seeded per DDI, so the membership rows
    // come out the same.
    //
    // It must stay *after* PRUNE_UNREPRESENTED_DDIS: ddi_split_membership.ddi_id
    // is ON DELETE CASCADE, so writing membership before the prune would lose
    // rows to it.
    //
    def external_methods = splitRows()
        .findAll { row -> !row.internal_test }
        .collectMany { row -> activeNegsets(row).collect { negset -> row.negsets[negset] } }

    external = BUILD_EXTERNAL_TEST(
        domainsplit_db,
        external_methods,
        'test',
    )

    //
    // Domain embeddings, published as one HDF5 per model rather than stored.
    //
    // Keyed on the *pruned* database, not the enriched one, so this whole branch
    // runs concurrently with enrichment. Same style of ordering invariant as
    // BUILD_EXTERNAL_TEST above: the export needs only the `(domain.id,
    // instance_id)` pairs of `domain_protein_map`, and nothing downstream of the
    // prune adds or removes a `domain` or a `domain_protein_map` row --
    // BUILD_EXTERNAL_TEST writes only `ddi_split_membership`, and with
    // INSERT_DOMAIN_PROTEIN_MAPPING deleted enrichment writes neither table.
    // SUBSET_SPLIT_DB copies `domain.id` verbatim and the prune deletes without
    // renumbering, so those pairs are exactly the ones in every published split
    // database of this run.
    //
    // The FASTA is ppi-splitting's `sequences.fasta` itself: it already holds the
    // cut domain sequences keyed by instance id, which is the key the export
    // needs, so nothing re-derives it.
    //
    generate_domain_embeddings(
        union_sequences,
        domainsplit_db,
    )

    //
    // Enrichment, on the master database, before any subsetting.
    //
    protein_domain_map = ingested.protein_domain_map

    ENRICH_DDI_DATABASE(
        external.domainsplit_db,
        input_pfam2go,
        input_uniprot_sequences,
        protein_domain_map,
        input_uniprot_go_terms,
        input_string,
        input_string_map,
    )

    enriched_db = ENRICH_DDI_DATABASE.out.domainsplit_db

    //
    // One database per (method, split). SUBSET_SPLIT_DB creates each output and
    // pulls the surviving rows out of the master read-only -- it does not clone
    // and delete, so a split's file costs only the rows it keeps.
    //
    def subset_targets = splitRows().collectMany { row ->
        activeNegsets(row).collectMany { negset ->
            def method = row.negsets[negset]
            def splits = row.internal_test ? splitsOf(row) : splitsOf(row) + ['test']
            splits.collect { split -> [method: method, split: split] }
        }
    }

    split_dbs = SUBSET_SPLIT_DB(
        channel.fromList(subset_targets).combine(enriched_db)
    )

    //
    // Collate and save software versions
    //
    // ppi-splitting captures no tool versions at all (no versions.yml, no
    // `topic:`), so there is nothing to mix in from it -- see the deviations
    // section of its integration plan.
    //
    ch_versions = ch_versions.mix(
        INIT_DOMAINSPLIT_DB.out.versions,
        COLLECT_DDI_DATA.out.versions,
        EXPORT_UNION_FAMILIES.out.versions,
        EXPORT_SPLIT_DDIS.out.versions,
        INGEST_SPLITS.out.versions,
        INSERT_EXTERNAL_SOURCES.out.versions,
        PRUNE_UNREPRESENTED_DDIS.out.versions,
        REPORT_DDI_ATTRITION.out.versions,
        generate_domain_embeddings.out.versions,
        ENRICH_DDI_DATABASE.out.versions,
        BUILD_EXTERNAL_TEST.out.versions,
        SUBSET_SPLIT_DB.out.versions,
    )

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'pipeline_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )

emit:
    // The master: enriched *and* carrying every split's membership rows,
    // including the external test set BUILD_EXTERNAL_TEST wrote upstream.
    domainsplit_db    = enriched_db
    split_db          = split_dbs.split_db
    domain_embeddings = generate_domain_embeddings.out.embeddings
    candidate_network = candidate_network
    source_conflicts  = inserted.conflicts
    ddi_attrition     = attrition.report
    bias_analysis     = ppi.multiqc_report
    dropped_families  = dropped_families
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
