/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SPLIT_COBINET_DATABASE -- split the CoBiNet DB into train/opt/test sets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Extracts protein/domain sequences, clusters with MMseqs2, then runs the
    three splitting strategies (random DDI, random denoise, minimal leakage)
    producing per-method/per-split SQLite databases.
----------------------------------------------------------------------------*/

include { RANDOM_DDI_SPLIT                                                                                          } from '../../../modules/local/random_ddi_split/main'
include { RANDOM_DENOISE_SPLIT                                                                                      } from '../../../modules/local/random_denoise_split/main'
include { RANDOM_DENOISE_SPLIT as RANDOM_DENOISE_SPLIT_2                                                            } from '../../../modules/local/random_denoise_split/main'
include { SPLIT_DATABASE                                                                                            } from '../../../modules/local/split_database/main'
include { EXTRACT_DOMAIN_SEQUENCES; EXTRACT_PROTEIN_SEQUENCES; MINIMAL_LEAKAGE_SPLIT_DOMAIN; MINIMAL_LEAKAGE_SPLIT_PROTEIN } from '../../../modules/local/minimal_leakage_split/main'
include { MMSEQS_EASYCLUSTER                                                                                        } from '../../../modules/nf-core/mmseqs/easycluster/main'


def map_split_files(split_id_files_ch, split_paths_ch, method) {
    return split_id_files_ch.flatMap()
        .join(
            split_paths_ch.flatten().map {
                [it.getName(), it]
            }
        )
        .map { [[id:("${method}_" + it[1]), split: it[1], method: method], it[2]] }
}

workflow SPLIT_COBINET_DATABASE {
    take:
    cobinet_db_ch

    main:
    // Extract protein and domain sequences for clustering
    protein_sequences = EXTRACT_PROTEIN_SEQUENCES(
        cobinet_db_ch
    ).protein_fasta
    domain_sequences = EXTRACT_DOMAIN_SEQUENCES(
        cobinet_db_ch
    ).domain_fasta

    // Cluster protein and domain sequences

    def cluster_input = protein_sequences.map { path ->
        tuple([id: 'protein'], path)
    }.concat(
        domain_sequences.map { path ->
            tuple([id: 'domain'], path)
        }
    )

    clusters = MMSEQS_EASYCLUSTER(cluster_input)

    def splits = [
        ["train", 0.6,],
        ["optimization", 0.2,],
        ["test", 0.2]
    ]

    (split_ddi_paths, split_ddi_id_files_ch) = RANDOM_DDI_SPLIT(
        cobinet_db_ch,
        Channel.of(splits)
    )

    (split_denoise_paths, split_denoise_id_files_ch) = RANDOM_DENOISE_SPLIT(
        cobinet_db_ch,
        0.1,
        false
    )

    (split_discovery_paths, split_discovery_id_files_ch) = RANDOM_DENOISE_SPLIT_2(
        cobinet_db_ch,
        0.1,
        true
    )

    (split_minimal_leakage_protein_paths, split_minimal_leakage_protein_id_files_ch) = MINIMAL_LEAKAGE_SPLIT_PROTEIN(
        cobinet_db_ch,
        splits,
        clusters.tsv.filter { it[0].id == "protein" }.map { it[1] }
    )

    (split_minimal_leakage_domain_paths, split_minimal_leakage_domain_id_files_ch) = MINIMAL_LEAKAGE_SPLIT_DOMAIN(
        cobinet_db_ch,
        splits,
        clusters.tsv.filter { it[0].id == "domain" }.map { it[1] }
    )

    split_ch = Channel.empty().mix(
        map_split_files(split_ddi_id_files_ch, split_ddi_paths, "random_ddi"),
        map_split_files(split_denoise_id_files_ch, split_denoise_paths, "random_denoise"),
        map_split_files(split_discovery_id_files_ch, split_discovery_paths, "random_discovery"),
        map_split_files(split_minimal_leakage_protein_id_files_ch, split_minimal_leakage_protein_paths, "minimal_leakage_protein"),
        map_split_files(split_minimal_leakage_domain_id_files_ch, split_minimal_leakage_domain_paths, "minimal_leakage_domain")
    )

    // Split the cobinet database based on the split DDI ID files
    split_db_ch = SPLIT_DATABASE(
        split_ch,
        cobinet_db_ch
    )


    emit:
    split_db = split_db_ch
}
