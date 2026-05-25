/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SPLIT_DOMAINSPLIT_DATABASE -- split the Domainsplit DB into train/opt/test
    sets using random and minimal-leakage strategies.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Extracts domain sequences, clusters with MMseqs2, then runs the
    splitting strategies (random DDI as biased baseline, spectral
    graph-partitioning minimal leakage on domains) producing per-split
    SQLite databases directly.
----------------------------------------------------------------------------*/

include { RANDOM_DDI_SPLIT                                                  } from '../../../modules/local/random_ddi_split/main'
include { EXTRACT_DOMAIN_SEQUENCES; MINIMAL_LEAKAGE_SPLIT_DOMAIN            } from '../../../modules/local/minimal_leakage_split/main'
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

    def splits = [
        ["train", 0.6],
        ["optimization", 0.2],
        ["test", 0.2]
    ]

    // Biased baseline: random DDI split (same proteins in train and test)
    RANDOM_DDI_SPLIT(
        domainsplit_db_ch,
        Channel.of(splits)
    )

    // Leakage-aware: spectral graph partitioning on domain clusters
    MINIMAL_LEAKAGE_SPLIT_DOMAIN(
        domainsplit_db_ch,
        splits,
        clusters.tsv.filter { it[0].id == "domain" }.map { it[1] }
    )

    split_ch = Channel.empty().mix(
        map_split_dbs(RANDOM_DDI_SPLIT.out.split_info, RANDOM_DDI_SPLIT.out.split_dbs, "random_ddi"),
        map_split_dbs(MINIMAL_LEAKAGE_SPLIT_DOMAIN.out.split_info, MINIMAL_LEAKAGE_SPLIT_DOMAIN.out.split_dbs, "minimal_leakage_domain")
    )

    emit:
    split_db = split_ch
}
