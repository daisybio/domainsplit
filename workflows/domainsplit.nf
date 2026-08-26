/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap            } from 'plugin/nf-schema'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_domainsplit_pipeline'
include { INIT_DOMAINSPLIT_DB         } from '../modules/local/init_domainsplit_db/main.nf'
include { COLLECT_DDI_DATA            } from '../subworkflows/local/collect_ddi_data/main.nf'
include { EXPORT_UNION_FAMILIES       } from '../modules/local/export_union_families/main.nf'
include { EXPORT_SPLIT_DDIS           } from '../modules/local/export_split_ddis/main.nf'
include { INGEST_SPLITS               } from '../subworkflows/local/ingest_splits/main.nf'
include { INSERT_EXTERNAL_SOURCES     } from '../modules/local/insert_external_sources/main.nf'
include { PRUNE_UNREPRESENTED_DDIS    } from '../modules/local/prune_unrepresented_ddis/main.nf'
include { generate_esm_embeddings     } from '../modules/local/esm_embeddings/main.nf'
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
    ppi-splitting's comma-separated `negative_sampling_method` fan-out; until
    that lands, `params.ppi_splitting_multi_negset = false` keeps one negative
    set per row and the `_hcni` directories are simply not produced.
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

// The negative sets a row actually produces on this run. Without the fan-out
// only the first one is asked for -- ppi-splitting would treat
// "ilp,ilp_candidates" as a single unknown method name and fall through to its
// degree-preserving default sampler, silently.
def activeNegsets(row) {
    return params.ppi_splitting_multi_negset ? row.negsets.keySet() as List : [row.negsets.keySet().first()]
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

    input_uniprot_id_mapping = file(params.url_uniprot_id_mapping)
    input_uniprot_go_terms   = file(params.url_uniprot_go_terms)
    input_uniprot_sequences  = file(params.url_uniprot_sequences)
    input_string             = file(params.url_string)
    input_pfam2go            = file(params.url_pfam2go)

    def prott5_file = file(params.url_uniprot_prott5_embeddings)

    empty_db = INIT_DOMAINSPLIT_DB().domainsplit_db

    //
    // 3did into the database; every other source parsed to a normalized TSV.
    //
    COLLECT_DDI_DATA(
        empty_db,
        params.url_3did,
        params.url_negatome,
        params.url_uniprot_swissprot_pfam,
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
    fasta_ch = params.pfam_fasta
        ? channel.value(file(params.pfam_fasta, checkIfExists: true))
        : channel.value([])
    cache_dir = params.interpro_cache ? file(params.interpro_cache).toAbsolutePath().toString() : ''
    if (cache_dir && !file(cache_dir).exists()) {
        log.warn "--interpro_cache ${cache_dir} does not exist yet; FETCH_DOMAIN_META will create it."
    }

    fetched = FETCH_DOMAIN_META(
        families.map { fams -> tuple([id: 'union'], fams) },
        clans_ch,
        fasta_ch,
        channel.value(cache_dir),
    )

    // The families Pfam had nothing usable for, one row per family with a
    // `reason` column. Taken from the channel, never read back from the
    // published path. Three reasons and they are not interchangeable:
    // `no_eligible_instances` is a consequence of `instance_tier`, which this
    // pipeline pins to `human_only`, so it fires in bulk and would bury the
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
    // actually resolved to instances. Under `instance_tier = human_only` a
    // family whose strata are all non-human ends with none, and its DDIs leave
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
    // The external test set: DDIs from PPIDM / single_domain_ppi / negatome,
    // none of which can have appeared in any split of any set. Written as the
    // `test` split of every external_test method directory.
    //
    // Deliberately *before* enrichment, and it is not an optimisation of
    // convenience: BUILD_EXTERNAL_TEST clones the master, and after enrichment
    // ~99% of that file is per-residue ProtT5/ESM blobs (590 MB of 594 MB at
    // `-profile test` scale, tens of GB in production) that this step never
    // reads. Running it here clones a few MB instead.
    //
    // Safe because the two write disjoint tables and the read sets do not move:
    // this step reads `domain_domain_interaction`, `domain` and
    // `domain_protein_map` WHERE instance_id IS NOT NULL, and writes only
    // `ddi_split_membership`. Enrichment never touches the first two or the
    // last; its one overlap is INSERT_DOMAIN_PROTEIN_MAPPING's upsert, which
    // leaves `instance_id` NULL on the rows it adds and updates only
    // domain_sequence / esm*_per_domain on the rows it hits -- so the set of
    // rows with a non-NULL instance_id, and their (domain_id, instance_id,
    // protein_id) triples, are identical either side of it. Sampling is seeded
    // per DDI, so the membership rows come out the same.
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
    // Enrichment, on the master database, before any subsetting.
    //
    protein_domain_map = ingested.protein_domain_map

    generate_esm_embeddings(
        input_uniprot_sequences,
        protein_domain_map,
    )

    ENRICH_DDI_DATABASE(
        external.domainsplit_db,
        input_pfam2go,
        input_uniprot_sequences,
        protein_domain_map,
        prott5_file,
        input_uniprot_go_terms,
        input_string,
        input_uniprot_id_mapping,
        generate_esm_embeddings.out.protein_embeddings,
        generate_esm_embeddings.out.domain_embeddings,
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
        generate_esm_embeddings.out.versions,
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
    candidate_network = candidate_network
    source_conflicts  = inserted.conflicts
    bias_analysis     = ppi.multiqc_report
    dropped_families  = dropped_families
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
