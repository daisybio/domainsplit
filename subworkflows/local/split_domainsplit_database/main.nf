/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SPLIT_DOMAINSPLIT_DATABASE -- split the Domainsplit DB into train/opt/test
    sets using random and minimal-leakage strategies.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Extracts domain sequences, clusters with MMseqs2, then runs the
    splitting strategies (random DDI as biased baseline, spectral
    graph-partitioning minimal leakage on domains) and materialises
    per-method/per-split SQLite databases.
----------------------------------------------------------------------------*/

include { RANDOM_DDI_SPLIT                                                  } from '../../../modules/local/random_ddi_split/main'
include { SPLIT_DATABASE                                                    } from '../../../modules/local/split_database/main'
include { EXTRACT_DOMAIN_SEQUENCES; MINIMAL_LEAKAGE_SPLIT_DOMAIN            } from '../../../modules/local/minimal_leakage_split/main'
include { MMSEQS_EASYCLUSTER                                                } from '../../../modules/nf-core/mmseqs/easycluster/main'


def map_split_files(split_id_files_ch, split_paths_ch, method) {
    return split_id_files_ch.flatMap()
        .join(
            split_paths_ch.flatten().map {
                [it.getName(), it]
            }
        )
        .map { [[id:("${method}_" + it[1]), split: it[1], method: method], it[2]] }
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
    (split_ddi_paths, split_ddi_id_files_ch) = RANDOM_DDI_SPLIT(
        domainsplit_db_ch,
        Channel.of(splits)
    )

    // Leakage-aware: spectral graph partitioning on domain clusters
    (split_minimal_leakage_domain_paths, split_minimal_leakage_domain_id_files_ch) = MINIMAL_LEAKAGE_SPLIT_DOMAIN(
        domainsplit_db_ch,
        splits,
        clusters.tsv.filter { it[0].id == "domain" }.map { it[1] }
    )

    split_ch = Channel.empty().mix(
        map_split_files(split_ddi_id_files_ch, split_ddi_paths, "random_ddi"),
        map_split_files(split_minimal_leakage_domain_id_files_ch, split_minimal_leakage_domain_paths, "minimal_leakage_domain")
    )

    // Materialise per-method/per-split SQLite databases
    SPLIT_DATABASE(
        split_ch,
        domainsplit_db_ch
    )

    emit:
    split_db = SPLIT_DATABASE.out.split_db
}
