process EXTRACT_PROTEIN_SEQUENCES {
    tag "proteins"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path "domainsplit.sqlite3"

    output:
    path "protein_sequences.fasta.gz", emit: protein_fasta
    path "versions.yml", emit: versions

    script:
    """
    set -o pipefail

    # x'0a' is a newline character in hexadecimal, which is used to separate the header and sequence in FASTA format
    sqlite3 domainsplit.sqlite3 "
            SELECT
                CONCAT('>', protein.id, x'0a', protein.sequence)
            FROM protein;
        " | \\
        gzip > \\
        protein_sequences.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqlite3: \$(sqlite3 --version | awk '{print \$1}')
    END_VERSIONS
    """
}

process EXTRACT_DOMAIN_SEQUENCES {
    tag "domains"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path "domainsplit.sqlite3"

    output:
    path "domain_sequences.fasta.gz", emit: domain_fasta
    path "versions.yml", emit: versions

    script:
    """
    set -o pipefail

    # x'0a' is a newline character in hexadecimal, which is used to separate the header and sequence in FASTA format
    sqlite3 domainsplit.sqlite3 "
            SELECT
                CONCAT('>', domain_id, '-', protein_id, x'0a', domain_sequence)
            FROM domain_protein_map;
        " | \\
        gzip > \\
        domain_sequences.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqlite3: \$(sqlite3 --version | awk '{print \$1}')
    END_VERSIONS
    """
}

// FUTURE WORK: MINIMAL_LEAKAGE_SPLIT_PROTEIN is not yet implemented.
// Body is intentionally `exit 1` so any subworkflow call hits a hard failure
// instead of producing silent empty splits. Mirror the domain-level annealing
// approach (see MINIMAL_LEAKAGE_SPLIT_DOMAIN below) using protein clusters from
// MMSEQS_EASYCLUSTER over uniprot sequences. Until implemented, callers must
// route around this process — the workflow already excludes it from the
// default split methods list.
process MINIMAL_LEAKAGE_SPLIT_PROTEIN {
    tag "minimal_leakage_protein"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path "domainsplit.sqlite3"
    val split_fractions  // e.g., [("train", 0.8), ("test", 0.2)]
    path ("protein_clusters.tsv")

    output:
    path output_files, emit: split_ddi_id_files
    val output_file_splits, emit: split_fractions

    script:
    def output_file_fraction_dict = [:]
    def output_file_splits = [:]

    split_fractions.each { name, fraction ->
        output_file_fraction_dict["${name}.txt"] = fraction
        output_file_splits["${name}.txt"] = name
    }

    def output_files = output_file_fraction_dict.keySet() //.collect { file_name -> file_name }

    def split_fraction_dict_str = output_file_fraction_dict.collect { k, v -> "'${k}': ${v}" }.join(", ")
    def split_fraction_dict_py = "{" + split_fraction_dict_str + "}"

    """
    echo "MINIMAL_LEAKAGE_SPLIT_PROTEIN is not yet implemented." >&2
    exit 1
    """
}

process MINIMAL_LEAKAGE_SPLIT_DOMAIN {
    tag "minimal_leakage_domain"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path "domainsplit.sqlite3"
    val split_fractions  // e.g., [("train", 0.8), ("test", 0.2)]
    path ("domain_clusters.tsv")

    output:
    path output_files, emit: split_ddi_id_files
    val output_file_splits, emit: split_fractions
    path "versions.yml", emit: versions

    script:
    def output_file_fraction_dict = [:]
    def output_file_splits = [:]

    split_fractions.each { name, fraction ->
        output_file_fraction_dict["${name}.txt"] = fraction
        output_file_splits["${name}.txt"] = name
    }

    def output_files = output_file_fraction_dict.keySet() //.collect { file_name -> file_name }

    def split_fraction_dict_str = output_file_fraction_dict.collect { k, v -> "'${k}': ${v}" }.join(", ")
    def split_fraction_dict_py = "{" + split_fraction_dict_str + "}"

    """
    #!/usr/bin/env python3

    import numpy as np
    import pandas as pd
    import random
    import sqlite3
    from functools import partial
    from tqdm import tqdm

    ANNEALING_STEPS = 3
    annealing_step_fractions = [(i + 1) / (ANNEALING_STEPS + 1) for i in range(ANNEALING_STEPS)]
    split_fractions = ${split_fraction_dict_py}
    print(f"Split fractions: {split_fractions}")

    # cluster sequences are defined as <domain_id>-<protein_id>
    cluster_df = pd.read_csv("domain_clusters.tsv", sep="\t", header=None, names=["centroid", "member"])
    # strip the protein_id from the sequence ids in the cluster dataframe
    cluster_df = cluster_df.map(lambda x: x.split("-")[0]).map(int)
    # Aggregate member domains for each centroid
    domain_clusters = cluster_df.groupby("centroid")["member"].apply(set).tolist()

    # Load DDI data from the database
    conn = sqlite3.connect("domainsplit.sqlite3")
    # Load the protein-domain mapping into a DataFrame (not used in the current implementation, but may be useful for future extensions)
    # pd_mapping = pd.read_sql('''
    #    SELECT domain_id, protein_id FROM domain_protein_map
    #    WHERE domain_id IS NOT NULL AND protein_id IS NOT NULL LIMIT 20;
    # ''', conn)

    ddi = pd.read_sql('''
        SELECT id AS ddi_id, domain_id_a, domain_id_b FROM domain_domain_interaction;
    ''', conn)

    conn.close()

    def get_new_ddis_for_domains(domains_in_split, new_domains):
        prefiltered_ddi = ddi[(ddi["domain_id_a"].isin(new_domains)) | (ddi["domain_id_b"].isin(new_domains))]
        new_ddis = {ddi_id
                    for ddi_id, domain_a, domain_b
                    in prefiltered_ddi.itertuples(index=False)
                    if (domain_a in domains_in_split or domain_b in domains_in_split)
                    }
        return new_ddis

    def get_split_ddis(split_name, new_domains=None):
        if new_domains is None:
            # compute from scratch (more expensive)
            return {ddi_id
                    for ddi_id, domain_a, domain_b
                    in ddi.itertuples(index=False)
                    if domain_a in split_domains[split_name] and domain_b in split_domains[split_name]
                    }
        else:
            # compute incrementally (more efficient)
            new_ddis = get_new_ddis_for_domains(split_domains[split_name], new_domains)
            return split_ddis[split_name].union(new_ddis)


    def cluster_rank_function(split_name, cluster_domains):
        # Calculate the total number of interactions that would be added to the split if this cluster were added
        # no need to normalize by split size since we are comparing clusters for the same split
        #return len(get_split_ddis(split_name, new_domains=cluster_domains))

        return len(get_new_ddis_for_domains(split_domains[split_name], cluster_domains))


    def split(list, fraction):
        split_index = int(len(list) * fraction)
        return list[:split_index], list[split_index:]

    split_domains = dict()
    split_ddis = dict()

    # as initialization, assign clusters to splits according to the annealing fraction, then iteratively assign remaining clusters to the split that would gain the most interactions (normalized by split size)
    for (split_name, fraction), clusters in zip(split_fractions.items(),
                                                np.array_split(domain_clusters, len(split_fractions))):
        split_domains[split_name] = set.union(*clusters)

    for annealing_step, annealing_fraction in enumerate(annealing_step_fractions):
        print(f"Starting annealing step {annealing_step + 1}/{ANNEALING_STEPS} with fraction {annealing_fraction:.2f}")
        print("Initial DDI counts: ")
        for split_name in split_fractions.keys():
            split_ddis[split_name] = get_split_ddis(split_name)
            print(f"\t{split_name}: {len(split_ddis[split_name])} interactions")

        random.shuffle(domain_clusters)
        clusters_to_keep, clusters_queue_list = split(domain_clusters, annealing_fraction)
        domains_to_keep = set.union(*clusters_to_keep)

        # Remove domains from splits that are not in the clusters to keep
        for split_name in split_fractions.keys():
            split_domains[split_name] = split_domains[split_name].intersection(domains_to_keep)
            split_ddis[split_name] = get_split_ddis(split_name)

        print("DDI counts after initialization:")
        for split_name in split_fractions.keys():
            print(f"\t{split_name}: {len(split_ddis[split_name])} interactions")

        # Keep clusters in a dict keyed by index so removal is O(1).
        # `clusters_queue.remove(cluster)` on a list was O(n) per iteration,
        # making the whole annealing loop O(n^2) for many clusters.
        clusters_queue = dict(enumerate(clusters_queue_list))
        tqdm_ = tqdm(total=len(domain_clusters), initial=len(clusters_to_keep), desc=f"Annealing step {annealing_step + 1}/{ANNEALING_STEPS}")
        while clusters_queue:
            # take the split that has the least number of interactions (normalized by split size)
            current_split_name = min(split_fractions.keys(),
                                     key=lambda split_name: len(split_ddis[split_name]) / split_fractions[split_name])

            # rank the clusters by the number of interactions they would add to the split
            cluster_idx, cluster = max(clusters_queue.items(),
                                       key=lambda kv: cluster_rank_function(current_split_name, kv[1]))
            del clusters_queue[cluster_idx]

            split_domains[current_split_name].update(cluster)
            split_ddis[current_split_name] = get_split_ddis(current_split_name, new_domains=cluster)

            tqdm_.update(1)

    print("Final DDI counts:")
    for output_file in split_fractions.keys():
        # output_file is e.g. "train.txt"; matches output_file_fraction_dict keys.
        with open(output_file, "w") as f:
            print(f"Writing split {output_file} with {len(split_domains[output_file])} domains and {len(split_ddis[output_file])} interactions to file...")
            f.write("ddi_id\\n")
            for domain_id in split_domains[output_file]:
                f.write(f"{domain_id}\\n")
            print("Done!")

    import sys as _sys
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {_sys.version.split()[0]}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
        f.write(f"    numpy: {np.__version__}\\n")
    """
}
