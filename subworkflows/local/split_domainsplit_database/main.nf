/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SPLIT_DOMAINSPLIT_DATABASE -- produce three split strategies:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      * random_ddi             biased baseline (random partition)
      * minimal_leakage_domain leakage-aware spectral partition
      * external_validation    leakage-aware train/validation on the "core"
                               sources (3did + high-conf non-PPI negatives),
                               plus an as-is test set from the held-out sources
                               (single-domain PPI, PPIDM, negatome).

    Domain sequences are extracted and clustered (MMseqs2) once; both the
    minimal-leakage and external-validation train/val partitions reuse the
    clusters.
----------------------------------------------------------------------------*/

include { RANDOM_DDI_SPLIT                                                  } from '../../../modules/local/random_ddi_split/main'
include { EXTRACT_DOMAIN_SEQUENCES; MINIMAL_LEAKAGE_SPLIT_DOMAIN            } from '../../../modules/local/minimal_leakage_split/main'
include { MINIMAL_LEAKAGE_SPLIT_DOMAIN as MINIMAL_LEAKAGE_SPLIT_TRAINVAL    } from '../../../modules/local/minimal_leakage_split/main'
include { SUBSET_DDIS_BY_SOURCE                                             } from '../../../modules/local/external_validation_split/main'
include { MMSEQS_EASYCLUSTER                                                } from '../../../modules/nf-core/mmseqs/easycluster/main'


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

    // The two methods that mimic a within-distribution evaluation use only the
    // "core" sources: 3did positives + high-confidence non-PPI negatives.
    // 'inferred_ppi_screen_negative' must stay in sync with the --source-label
    // default in bin/insert_ppi_negative_selection.py.
    def core_sources = ['3did', 'inferred_ppi_screen_negative']

    // External-validation test set: held-out sources placed as is.
    def test_sources = [
        'single_domain_ppi',
        'PPIDM_Bronze', 'PPIDM_Silver', 'PPIDM_Gold',
        'negatome',
    ]

    // Biased baseline: random DDI split (same proteins in train and test)
    RANDOM_DDI_SPLIT(
        domainsplit_db_ch,
        Channel.of(splits),
        core_sources
    )

    // Leakage-aware: spectral graph partitioning on domain clusters
    MINIMAL_LEAKAGE_SPLIT_DOMAIN(
        domainsplit_db_ch,
        splits,
        clusters_tsv,
        core_sources
    )

    // External validation: leakage-free train/validation on core sources ...
    def trainval_splits = [
        ["train", 0.8],
        ["validation", 0.2]
    ]
    MINIMAL_LEAKAGE_SPLIT_TRAINVAL(
        domainsplit_db_ch,
        trainval_splits,
        clusters_tsv,
        core_sources
    )

    // ... plus an as-is test set from the held-out sources
    SUBSET_DDIS_BY_SOURCE(
        domainsplit_db_ch,
        test_sources,
        "test"
    )

    split_ch = Channel.empty().mix(
        map_split_dbs(RANDOM_DDI_SPLIT.out.split_info, RANDOM_DDI_SPLIT.out.split_dbs, "random_ddi"),
        map_split_dbs(MINIMAL_LEAKAGE_SPLIT_DOMAIN.out.split_info, MINIMAL_LEAKAGE_SPLIT_DOMAIN.out.split_dbs, "minimal_leakage_domain"),
        map_split_dbs(MINIMAL_LEAKAGE_SPLIT_TRAINVAL.out.split_info, MINIMAL_LEAKAGE_SPLIT_TRAINVAL.out.split_dbs, "external_validation"),
        map_split_dbs(SUBSET_DDIS_BY_SOURCE.out.split_info, SUBSET_DDIS_BY_SOURCE.out.split_dbs, "external_validation")
    )

    emit:
    split_db = split_ch
}
