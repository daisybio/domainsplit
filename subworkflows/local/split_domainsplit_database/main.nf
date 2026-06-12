/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SPLIT_DOMAINSPLIT_DATABASE -- produce the split strategies, each run ONCE
    PER NEGATIVE-DDI METHOD ("deletion" and "random_addition") so the two
    methods' core sources stay isolated.  Strategies:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      * random_ddi             biased baseline (random partition)
      * minimal_leakage_domain leakage-aware spectral partition
      * external_validation    leakage-aware train/validation on the "core"
                               sources (3did copy + the method's PPI-screen
                               negatives), plus an as-is test set from the
                               held-out sources (single-domain PPI, PPIDM,
                               negatome).  The held-out test set is
                               method-independent, so it is built once and routed
                               into both external_validation_* folders.

    Each strategy therefore yields two method folders, e.g.
    random_ddi_deletion / random_ddi_random_addition -> 6 method folders total.

    Domain sequences are extracted and clustered (MMseqs2) once; every
    leakage-aware partition reuses the clusters.
----------------------------------------------------------------------------*/

include { RANDOM_DDI_SPLIT as RANDOM_DDI_SPLIT_DEL          } from '../../../modules/local/random_ddi_split/main'
include { RANDOM_DDI_SPLIT as RANDOM_DDI_SPLIT_RAND         } from '../../../modules/local/random_ddi_split/main'
include { EXTRACT_DOMAIN_SEQUENCES                          } from '../../../modules/local/minimal_leakage_split/main'
include { MINIMAL_LEAKAGE_SPLIT_DOMAIN as MLS_DOMAIN_DEL    } from '../../../modules/local/minimal_leakage_split/main'
include { MINIMAL_LEAKAGE_SPLIT_DOMAIN as MLS_DOMAIN_RAND   } from '../../../modules/local/minimal_leakage_split/main'
include { MINIMAL_LEAKAGE_SPLIT_DOMAIN as MLS_TRAINVAL_DEL  } from '../../../modules/local/minimal_leakage_split/main'
include { MINIMAL_LEAKAGE_SPLIT_DOMAIN as MLS_TRAINVAL_RAND } from '../../../modules/local/minimal_leakage_split/main'
include { SUBSET_DDIS_BY_SOURCE                             } from '../../../modules/local/external_validation_split/main'
include { MMSEQS_EASYCLUSTER                                } from '../../../modules/nf-core/mmseqs/easycluster/main'


def map_split_dbs(split_info_ch, split_dbs_ch, method) {
    return split_dbs_ch.flatten().map { db ->
        [db.getName(), db]
    }.join(
        split_info_ch.flatMap()
    ).map { [[ id: "${method}_${it[2]}", split: it[2], method: method ], it[1]] }
}

workflow SPLIT_DOMAINSPLIT_DATABASE {
    take:
    domainsplit_db_ch

    main:
    // Extract domain sequences for clustering
    domain_sequences = EXTRACT_DOMAIN_SEQUENCES(
        domainsplit_db_ch
    ).domain_fasta

    // Cluster domain sequences (groups sequence-similar domains)
    def cluster_input = domain_sequences.map { path ->
        tuple([id: 'domain'], path)
    }

    clusters = MMSEQS_EASYCLUSTER(cluster_input)
    def clusters_tsv = clusters.tsv.filter { it[0].id == "domain" }.map { it[1] }

    def splits = [
        ["train", 0.6],
        ["optimization", 0.2],
        ["test", 0.2]
    ]
    def trainval_splits = [
        ["train", 0.8],
        ["validation", 0.2]
    ]

    // Core sources per negative-DDI method: the method's 3did positive copy plus
    // its PPI-screen negatives.  Must stay in sync with the source labels written
    // by bin/insert_ppi_negative_selection.py.
    def core_deletion = [
        '3did_deletion',
        'inferred_ppi_screen_negative_for_deletion',
    ]
    def core_random_addition = [
        '3did_random_addition',
        'inferred_ppi_screen_negative_for_random_addition',
    ]

    // External-validation test set: held-out sources placed as is
    // (method-independent).
    def test_sources = [
        'single_domain_ppi',
        'PPIDM_Bronze', 'PPIDM_Silver', 'PPIDM_Gold',
        'negatome',
    ]

    // Biased baseline: random DDI split (same proteins in train and test)
    RANDOM_DDI_SPLIT_DEL(domainsplit_db_ch, splits, core_deletion)
    RANDOM_DDI_SPLIT_RAND(domainsplit_db_ch, splits, core_random_addition)

    // Leakage-aware: spectral graph partitioning on domain clusters
    MLS_DOMAIN_DEL(domainsplit_db_ch, splits, clusters_tsv, core_deletion)
    MLS_DOMAIN_RAND(domainsplit_db_ch, splits, clusters_tsv, core_random_addition)

    // External validation: leakage-free train/validation on core sources ...
    MLS_TRAINVAL_DEL(domainsplit_db_ch, trainval_splits, clusters_tsv, core_deletion)
    MLS_TRAINVAL_RAND(domainsplit_db_ch, trainval_splits, clusters_tsv, core_random_addition)

    // ... plus an as-is test set from the held-out sources (shared by both
    //     external_validation_* methods).
    SUBSET_DDIS_BY_SOURCE(
        domainsplit_db_ch,
        test_sources,
        "test"
    )

    split_ch = Channel.empty().mix(
        map_split_dbs(RANDOM_DDI_SPLIT_DEL.out.split_info,  RANDOM_DDI_SPLIT_DEL.out.split_dbs,  "random_ddi_deletion"),
        map_split_dbs(RANDOM_DDI_SPLIT_RAND.out.split_info, RANDOM_DDI_SPLIT_RAND.out.split_dbs, "random_ddi_random_addition"),
        map_split_dbs(MLS_DOMAIN_DEL.out.split_info,  MLS_DOMAIN_DEL.out.split_dbs,  "minimal_leakage_domain_deletion"),
        map_split_dbs(MLS_DOMAIN_RAND.out.split_info, MLS_DOMAIN_RAND.out.split_dbs, "minimal_leakage_domain_random_addition"),
        map_split_dbs(MLS_TRAINVAL_DEL.out.split_info,  MLS_TRAINVAL_DEL.out.split_dbs,  "external_validation_deletion"),
        map_split_dbs(MLS_TRAINVAL_RAND.out.split_info, MLS_TRAINVAL_RAND.out.split_dbs, "external_validation_random_addition"),
        map_split_dbs(SUBSET_DDIS_BY_SOURCE.out.split_info, SUBSET_DDIS_BY_SOURCE.out.split_dbs, "external_validation_deletion"),
        map_split_dbs(SUBSET_DDIS_BY_SOURCE.out.split_info, SUBSET_DDIS_BY_SOURCE.out.split_dbs, "external_validation_random_addition")
    )

    // NB: MMSEQS_EASYCLUSTER (nf-core) reports its version via the `versions`
    // channel topic, not an `emit: versions` output, so it is not mixed here.
    ch_versions = Channel.empty().mix(
        EXTRACT_DOMAIN_SEQUENCES.out.versions,
        RANDOM_DDI_SPLIT_DEL.out.versions,
        RANDOM_DDI_SPLIT_RAND.out.versions,
        MLS_DOMAIN_DEL.out.versions,
        MLS_DOMAIN_RAND.out.versions,
        MLS_TRAINVAL_DEL.out.versions,
        MLS_TRAINVAL_RAND.out.versions,
        SUBSET_DDIS_BY_SOURCE.out.versions,
    )

    emit:
    split_db = split_ch
    versions = ch_versions
}
