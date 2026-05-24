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
    output_file_fraction_dict = [:]
    output_file_splits = [:]

    split_fractions.each { name, fraction ->
        output_file_fraction_dict["${name}.txt"] = fraction
        output_file_splits["${name}.txt"] = name
    }

    output_files = output_file_fraction_dict.keySet() as List

    def split_fraction_dict_str = output_file_fraction_dict.collect { k, v -> "'${k}': ${v}" }.join(", ")
    def split_fraction_dict_py = "{" + split_fraction_dict_str + "}"

    """
    #!/usr/bin/env python3
    \"\"\"
    Minimal-leakage domain-level splitting via spectral graph partitioning.

    Algorithm (adapted from KaHIP/KaFFPa as used in Bernett et al. 2024,
    "Cracking the black box of deep sequence-based PPI prediction"):

    1. Build a weighted cluster-DDI graph where nodes = MMseqs2 domain
       clusters and edge weights = number of DDIs connecting domains
       across two clusters.
    2. Compute the graph Laplacian and its smallest non-trivial
       eigenvectors (Fiedler vectors) for a spectral embedding.
    3. Run weighted k-means in spectral space to get a balanced initial
       k-way partition respecting target split fractions.
    4. Refine with Kernighan-Lin local search: greedily move clusters
       between partitions when the move reduces the total cut (DDIs
       lost) while maintaining balance.
    5. Output DDI IDs per split — a DDI is assigned to a split iff both
       its domains land in the same partition.
    \"\"\"

    import numpy as np
    import sqlite3
    from collections import defaultdict

    split_fractions = ${split_fraction_dict_py}
    print(f"Split fractions: {split_fractions}")
    split_names = list(split_fractions.keys())
    k = len(split_names)

    # ── Load cluster data ────────────────────────────────────────────
    # Cluster TSV format: centroid<TAB>member, IDs are "domain_id-protein_id"
    domain_to_cluster = {}
    cluster_domains = defaultdict(set)
    with open("domain_clusters.tsv") as fh:
        for line in fh:
            parts = line.strip().split("\\t")
            if len(parts) != 2:
                continue
            centroid_domain = int(parts[0].split("-")[0])
            member_domain = int(parts[1].split("-")[0])
            cluster_domains[centroid_domain].add(member_domain)
            domain_to_cluster[member_domain] = centroid_domain

    # ── Load DDI data ────────────────────────────────────────────────
    conn = sqlite3.connect("domainsplit.sqlite3")
    ddi_rows = conn.execute(
        "SELECT id, domain_id_a, domain_id_b FROM domain_domain_interaction"
    ).fetchall()
    conn.close()
    print(f"Loaded {len(ddi_rows)} DDIs")

    # Domains not in any cluster → singleton clusters
    unclustered = 0
    for _, da, db in ddi_rows:
        for d in (da, db):
            if d not in domain_to_cluster:
                domain_to_cluster[d] = d
                cluster_domains[d] = {d}
                unclustered += 1
    if unclustered:
        print(f"Created {unclustered} singleton clusters for unclustered domains")

    cluster_ids = sorted(cluster_domains.keys())
    cluster_index = {cid: i for i, cid in enumerate(cluster_ids)}
    n_clusters = len(cluster_ids)
    print(f"Total clusters: {n_clusters}")

    # ── Build cluster-level DDI graph ────────────────────────────────
    cluster_pair_ddis = defaultdict(set)
    for ddi_id, da, db in ddi_rows:
        ca = domain_to_cluster.get(da)
        cb = domain_to_cluster.get(db)
        if ca is None or cb is None:
            continue
        key = (min(ca, cb), max(ca, cb))
        cluster_pair_ddis[key].add(ddi_id)

    print(f"Cluster graph: {n_clusters} nodes, {len(cluster_pair_ddis)} edges")

    cluster_weights = np.array(
        [len(cluster_domains[cid]) for cid in cluster_ids], dtype=float
    )
    total_weight = cluster_weights.sum()
    target_weights = np.array([split_fractions[s] for s in split_names]) * total_weight

    # ── Step 1: Spectral initialization ──────────────────────────────
    W = np.zeros((n_clusters, n_clusters))
    for (ci, cj), ddis in cluster_pair_ddis.items():
        i, j = cluster_index[ci], cluster_index[cj]
        w = len(ddis)
        if i == j:
            W[i, j] = w
        else:
            W[i, j] = w
            W[j, i] = w

    deg = W.sum(axis=1)
    L = np.diag(deg) - W

    print(f"Computing Laplacian eigenvectors ({n_clusters}x{n_clusters})...")
    eigenvalues, eigenvectors = np.linalg.eigh(L)
    spectral_coords = eigenvectors[:, 1:k]

    # ── Weighted k-means in spectral space ───────────────────────────
    def weighted_kmeans(coords, weights, k, target_wts, max_iter=100, seed=42):
        rng = np.random.RandomState(seed)
        n = len(coords)
        probs = weights / weights.sum()
        centers = coords[rng.choice(n, size=k, replace=False, p=probs)].copy()
        labels = np.zeros(n, dtype=int)

        for _ in range(max_iter):
            dists = np.array(
                [np.linalg.norm(coords - c, axis=1) for c in centers]
            ).T
            new_labels = np.argmin(dists, axis=1)

            # Greedy rebalancing
            part_wts = np.array(
                [weights[new_labels == j].sum() for j in range(k)]
            )
            for _ in range(20):
                excess = part_wts - target_wts
                over = int(np.argmax(excess))
                under = int(np.argmin(excess))
                if excess[over] <= target_wts[over] * 0.1:
                    break
                mask = new_labels == over
                if not mask.any():
                    break
                d2u = np.linalg.norm(coords[mask] - centers[under], axis=1)
                for idx in np.where(mask)[0][np.argsort(d2u)]:
                    new_labels[idx] = under
                    part_wts[over] -= weights[idx]
                    part_wts[under] += weights[idx]
                    if part_wts[over] <= target_wts[over] * 1.1:
                        break

            if np.array_equal(labels, new_labels):
                break
            labels = new_labels
            for j in range(k):
                m = labels == j
                if m.any():
                    centers[j] = np.average(coords[m], weights=weights[m], axis=0)

        return labels

    print("Weighted k-means in spectral space...")
    labels = weighted_kmeans(spectral_coords, cluster_weights, k, target_weights)

    partition = {name: set() for name in split_names}
    for i, cid in enumerate(cluster_ids):
        partition[split_names[labels[i]]].add(cid)

    # ── Step 2: Kernighan-Lin local search ───────────────────────────
    def compute_cut(partition):
        c2p = {}
        for name, clusters in partition.items():
            for c in clusters:
                c2p[c] = name
        cut = 0
        internal = 0
        for (ci, cj), ddis in cluster_pair_ddis.items():
            w = len(ddis)
            if ci == cj:
                internal += w
            elif c2p.get(ci) == c2p.get(cj):
                internal += w
            else:
                cut += w
        return cut, internal

    cut, internal = compute_cut(partition)
    print(f"Spectral init: {internal} DDIs preserved, {cut} cut")

    print("Kernighan-Lin refinement...")
    tolerance = 0.05
    max_kl_iter = 50

    for iteration in range(max_kl_iter):
        c2p = {}
        for name, clusters in partition.items():
            for c in clusters:
                c2p[c] = name

        cpw = defaultdict(lambda: defaultdict(int))
        for (ci, cj), ddis in cluster_pair_ddis.items():
            if ci == cj:
                continue
            w = len(ddis)
            cpw[ci][c2p.get(cj, "")] += w
            cpw[cj][c2p.get(ci, "")] += w

        part_wts = {
            name: sum(len(cluster_domains[c]) for c in clusters)
            for name, clusters in partition.items()
        }

        best_gain = 0
        best_move = None

        for cid in cluster_ids:
            src = c2p.get(cid)
            if src is None:
                continue
            cw = len(cluster_domains[cid])

            for dst in split_names:
                if dst == src:
                    continue
                new_src_wt = part_wts[src] - cw
                new_dst_wt = part_wts[dst] + cw
                if (new_src_wt < split_fractions[src] * total_weight * (1 - tolerance) or
                        new_dst_wt > split_fractions[dst] * total_weight * (1 + tolerance)):
                    continue
                gain = cpw[cid].get(dst, 0) - cpw[cid].get(src, 0)
                if gain > best_gain:
                    best_gain = gain
                    best_move = (cid, src, dst)

        if not best_move:
            print(f"  Converged at iteration {iteration + 1}")
            break

        cid, src, dst = best_move
        partition[src].remove(cid)
        partition[dst].add(cid)

        if (iteration + 1) % 10 == 0:
            c, i = compute_cut(partition)
            print(f"  Iteration {iteration + 1}: {i} preserved, {c} cut")

    cut, internal = compute_cut(partition)
    print(f"Final: {internal} DDIs preserved, {cut} cut")

    # ── Step 3: Assign DDIs to splits ────────────────────────────────
    c2p = {}
    for name, clusters in partition.items():
        for c in clusters:
            c2p[c] = name

    split_ddis = {name: [] for name in split_names}
    skipped = 0
    for ddi_id, da, db in ddi_rows:
        ca = domain_to_cluster.get(da)
        cb = domain_to_cluster.get(db)
        if ca is None or cb is None:
            skipped += 1
            continue
        pa, pb = c2p.get(ca), c2p.get(cb)
        if pa == pb and pa is not None:
            split_ddis[pa].append(ddi_id)

    total_assigned = sum(len(v) for v in split_ddis.values())
    print(f"Assigned {total_assigned}/{len(ddi_rows)} DDIs "
          f"({100 * total_assigned / max(len(ddi_rows), 1):.1f}%)")

    # ── Write output ─────────────────────────────────────────────────
    for output_file in split_fractions.keys():
        with open(output_file, "w") as f:
            f.write("ddi_id\\n")
            for ddi_id in sorted(split_ddis[output_file]):
                f.write(f"{ddi_id}\\n")
        n_domains = sum(len(cluster_domains[c]) for c in partition[output_file])
        frac = split_fractions[output_file]
        print(f"  {output_file}: {len(split_ddis[output_file])} DDIs, "
              f"{n_domains} domains "
              f"({100 * n_domains / total_weight:.1f}%, target {100 * frac:.0f}%)")

    import sys as _sys
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {_sys.version.split()[0]}\\n")
        f.write(f"    numpy: {np.__version__}\\n")
    """
}
