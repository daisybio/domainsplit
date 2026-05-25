#!/usr/bin/env python3
"""
Build negative DDIs from a Y2H/MS PPI parquet and append them to the
domainsplit SQLite, restricted to Pfam domains already present in positive
DDIs, ranked by how often each Pfam-pair co-occurs across the filtered PPI
set, capped so the new-source count equals (n_positive - n_negatome).
"""

import argparse
import gzip
import itertools
import json
import os
import sqlite3
import sys
from collections import defaultdict

import math

import numpy as np
import pyarrow.parquet as pq


TAG = "[ppi_neg]"
BATCH_SIZE = 500_000
REQUIRED_COLUMNS = ["gene_name_bait", "gene_name_prey", "n_tested"]


def log(msg):
    print(f"{TAG} {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True)
    p.add_argument("--parquet", required=True)
    p.add_argument("--idmapping", required=True, help="UniProt idmapping .dat or .dat.gz")
    p.add_argument("--pfam-stockholm", required=True,
                   help="Pfam-A.full.gz Stockholm alignment file")
    p.add_argument("--pfam-mapping-out", required=True,
                   help="Output path for UniProt→Pfam JSON mapping")
    p.add_argument("--min-n-tested", type=int, required=True)
    p.add_argument("--source-label", required=True)
    p.add_argument(
        "--sampling-strategy",
        choices=["frequency", "degree_matched"],
        default="degree_matched",
        help="'frequency' = top-N by co-occurrence (old behavior). "
             "'degree_matched' = sample to match positive degree distribution.",
    )
    return p.parse_args()


def open_idmapping(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_pfam_stockholm(path):
    if not os.path.isfile(path):
        sys.exit(f"{TAG} Stockholm file not found: {path}")

    uniprot_to_pfams = defaultdict(set)
    current_pfam = None
    n_entries = 0

    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#=GF AC"):
                current_pfam = line.split()[-1].split(".")[0]
                n_entries += 1
            elif line.startswith("#=GS") and "\tAC\t" not in line and " AC " in line:
                parts = line.split()
                try:
                    ac_idx = parts.index("AC")
                    uniprot = parts[ac_idx + 1].split(".")[0]
                    if current_pfam:
                        uniprot_to_pfams[uniprot].add(current_pfam)
                except (ValueError, IndexError):
                    continue
            elif line.startswith("//"):
                current_pfam = None

    log(f"parsed {n_entries} Pfam entries, {len(uniprot_to_pfams)} UniProt mappings")
    return uniprot_to_pfams


def load_positive_pfams(conn):
    cur = conn.execute(
        "SELECT DISTINCT d.pfam_id "
        "FROM domain AS d JOIN domain_domain_interaction AS ddi "
        "  ON d.id IN (ddi.domain_id_a, ddi.domain_id_b) "
        "WHERE ddi.negative = 0"
    )
    return {row[0] for row in cur}


def load_existing_pairs(conn):
    cur = conn.execute(
        "SELECT da.pfam_id, db.pfam_id "
        "FROM domain_domain_interaction AS ddi "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b"
    )
    return {tuple(sorted((a, b))) for a, b in cur}


def stream_idmapping(path, gene_set):
    gene_to_uniprots = defaultdict(set)
    with open_idmapping(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            uniprot, kind, value = parts
            if kind == "Gene_Name" or kind == "Gene_Synonym":
                if value in gene_set:
                    gene_to_uniprots[value].add(uniprot)
    return gene_to_uniprots


def _validate_columns(parquet_schema):
    available = set(parquet_schema.names)
    missing = set(REQUIRED_COLUMNS) - available
    if missing:
        sys.exit(f"{TAG} parquet missing required columns: {missing}")


def _collect_genes_and_pairs(parquet_path, min_n_tested):
    """Stream parquet in batches. Returns (unique_genes, filtered bait/prey lists)."""
    unique_genes = set()
    baits = []
    preys = []
    n_input = 0

    for batch in pq.ParquetFile(parquet_path).iter_batches(
        batch_size=BATCH_SIZE, columns=REQUIRED_COLUMNS
    ):
        tbl = batch.to_pydict()
        n_input += len(tbl["n_tested"])

        for bait, prey, n_tested in zip(
            tbl["gene_name_bait"], tbl["gene_name_prey"], tbl["n_tested"]
        ):
            if n_tested is None or n_tested < min_n_tested:
                continue
            if bait is None or prey is None:
                continue
            unique_genes.add(bait)
            unique_genes.add(prey)
            baits.append(bait)
            preys.append(prey)

    return n_input, unique_genes, baits, preys


def _compute_positive_degree(conn):
    """Per-domain degree in the positive DDI set."""
    rows = conn.execute(
        "SELECT domain_id_a, domain_id_b FROM domain_domain_interaction WHERE negative = 0"
    ).fetchall()
    deg = defaultdict(int)
    for a, b in rows:
        deg[a] += 1
        deg[b] += 1
    return deg


def select_degree_matched(fresh_candidates, pos_degree, n_take):
    """Select negatives matching the positive degree distribution shape.

    Bins positive degrees into log-spaced histogram buckets, then greedily
    picks candidates whose domains fall in under-represented bins, weighted
    by PPI co-occurrence count to preserve biological signal.
    """
    if not fresh_candidates or n_take <= 0:
        return []

    pos_vals = np.array(list(pos_degree.values()), dtype=float)
    if len(pos_vals) == 0:
        return fresh_candidates[:n_take]

    n_bins = min(10, max(3, int(np.sqrt(len(pos_vals)))))
    bin_edges = np.logspace(0, np.log10(max(pos_vals.max(), 1) + 1), n_bins + 1)
    bin_edges[0] = 0

    pos_hist, _ = np.histogram(pos_vals, bins=bin_edges)
    total_pos_pairs = sum(pos_hist)
    if total_pos_pairs == 0:
        return fresh_candidates[:n_take]

    target_per_bin = np.array([
        max(1, int(round(n_take * count / total_pos_pairs)))
        for count in pos_hist
    ], dtype=float)
    # Scale so sum matches n_take
    target_per_bin = target_per_bin * (n_take / target_per_bin.sum())

    filled_per_bin = np.zeros(n_bins)

    def domain_bin(pfam_id):
        deg = pos_degree.get(pfam_id, 0)
        idx = int(np.searchsorted(bin_edges, deg, side="right")) - 1
        return max(0, min(n_bins - 1, idx))

    def pair_score(pfam_a, pfam_b, count):
        bin_a = domain_bin(pfam_a)
        bin_b = domain_bin(pfam_b)
        deficit_a = max(0, target_per_bin[bin_a] - filled_per_bin[bin_a])
        deficit_b = max(0, target_per_bin[bin_b] - filled_per_bin[bin_b])
        return (deficit_a + deficit_b) * math.log1p(count)

    scored = []
    for (pfam_a, pfam_b), count in fresh_candidates:
        s = pair_score(pfam_a, pfam_b, count)
        scored.append(((pfam_a, pfam_b), count, s))

    scored.sort(key=lambda x: x[2], reverse=True)

    chosen = []
    for (pfam_a, pfam_b), count, _ in scored:
        if len(chosen) >= n_take:
            break
        chosen.append(((pfam_a, pfam_b), count))
        filled_per_bin[domain_bin(pfam_a)] += 1
        filled_per_bin[domain_bin(pfam_b)] += 1

    log(f"degree_matched: target_per_bin = {target_per_bin.astype(int).tolist()}")
    log(f"degree_matched: filled_per_bin = {filled_per_bin.astype(int).tolist()}")
    return chosen


def main():
    args = parse_args()

    pf = pq.ParquetFile(args.parquet)
    log(f"parquet {args.parquet}: {pf.metadata.num_rows} rows, "
        f"{pf.metadata.num_columns} cols, streaming in batches of {BATCH_SIZE}")
    _validate_columns(pf.schema_arrow)

    log("pass 1: collecting genes from parquet (batched)")
    n_input, unique_genes, baits, preys = _collect_genes_and_pairs(
        args.parquet, args.min_n_tested
    )
    n_after = len(baits)
    log(f"n_input_ppis = {n_input}")
    log(f"n_ppis_after_n_tested_filter (>= {args.min_n_tested}) = {n_after}")
    log(f"n_unique_genes = {len(unique_genes)}")

    log(f"parsing Stockholm file {args.pfam_stockholm}")
    uniprot_to_pfams = parse_pfam_stockholm(args.pfam_stockholm)

    log(f"writing UniProt→Pfam mapping to {args.pfam_mapping_out}")
    with open(args.pfam_mapping_out, "w") as fh:
        json.dump(
            {k: sorted(v) for k, v in uniprot_to_pfams.items()},
            fh,
        )

    log(f"streaming idmapping {args.idmapping}")
    gene_to_uniprots = stream_idmapping(args.idmapping, unique_genes)
    log(f"n_mapped_genes = {len(gene_to_uniprots)}")
    log(f"n_mapped_uniprots = {sum(len(v) for v in gene_to_uniprots.values())}")
    n_pfam_unique = len({p for s in uniprot_to_pfams.values() for p in s})
    log(f"n_pfam_domains_for_input_proteins = {n_pfam_unique}")

    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    pos_pfam = load_positive_pfams(conn)
    log(f"n_positive_pfams = {len(pos_pfam)}")

    existing_pairs = load_existing_pairs(conn)
    log(f"n_existing_ddis = {len(existing_pairs)}")

    def row_pfams(gene):
        return {
            pfam
            for uniprot in gene_to_uniprots.get(gene, ())
            for pfam in uniprot_to_pfams.get(uniprot, ())
        }

    candidate_counts = defaultdict(int)
    n_rows_with_pairs = 0
    for bait, prey in zip(baits, preys):
        bait_pfams = row_pfams(bait) & pos_pfam
        prey_pfams = row_pfams(prey) & pos_pfam
        if not bait_pfams or not prey_pfams:
            continue
        row_pairs = set()
        for a, b in itertools.product(bait_pfams, prey_pfams):
            if a == b:
                continue
            row_pairs.add(tuple(sorted((a, b))))
        if not row_pairs:
            continue
        n_rows_with_pairs += 1
        for key in row_pairs:
            candidate_counts[key] += 1

    del baits, preys

    log(f"n_rows_yielding_pairs = {n_rows_with_pairs}")
    log(f"n_unique_pfam_ddi_candidates = {len(candidate_counts)}")
    if candidate_counts:
        most_common_key, most_common_count = max(
            candidate_counts.items(), key=lambda kv: kv[1]
        )
        log(
            f"most_common_ddi = {most_common_key[0]}-{most_common_key[1]} "
            f"(observed in {most_common_count} PPI rows)"
        )

    n_positive = conn.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction WHERE negative = 0"
    ).fetchone()[0]
    n_negatome = conn.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction "
        "WHERE negative = 1 AND source = 'negatome'"
    ).fetchone()[0]
    log(f"n_positive_ddis_in_db = {n_positive}")
    log(f"n_negatome_negatives_in_db = {n_negatome}")
    n_take = max(0, n_positive - n_negatome)
    log(f"n_take (target for source='{args.source_label}') = {n_take}")

    fresh_candidates = [
        (key, count)
        for key, count in candidate_counts.items()
        if key not in existing_pairs
    ]
    fresh_candidates.sort(key=lambda kv: kv[1], reverse=True)
    log(f"n_fresh_candidates_after_dedup = {len(fresh_candidates)}")

    if args.sampling_strategy == "degree_matched":
        log("using degree-matched sampling strategy")
        pos_degree = _compute_positive_degree(conn)
        chosen = select_degree_matched(fresh_candidates, pos_degree, n_take)
    else:
        log("using frequency-ranked sampling strategy")
        chosen = fresh_candidates[:n_take]
    log(f"n_chosen = {len(chosen)}")

    if chosen:
        # Pre-load pfam_id → domain.id mapping to avoid per-row subqueries
        pfam_to_domain_ids = defaultdict(list)
        for did, pfam in conn.execute("SELECT id, pfam_id FROM domain"):
            pfam_to_domain_ids[pfam].append(did)
        log(f"loaded {len(pfam_to_domain_ids)} pfam→domain mappings")

        insert_rows = []
        for (pfam_a, pfam_b), _ in chosen:
            for d_a in pfam_to_domain_ids.get(pfam_a, ()):
                for d_b in pfam_to_domain_ids.get(pfam_b, ()):
                    insert_rows.append((d_a, d_b, True, args.source_label))

        conn.executemany(
            "INSERT OR IGNORE INTO domain_domain_interaction"
            "(domain_id_a, domain_id_b, negative, source) "
            "VALUES (?, ?, ?, ?)",
            insert_rows,
        )
        conn.commit()
        log(f"batch-inserted {len(insert_rows)} rows")

    n_inserted = conn.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction WHERE source = ?",
        (args.source_label,),
    ).fetchone()[0]
    log(f"n_inserted_for_source = {n_inserted}")
    conn.close()


if __name__ == "__main__":
    main()
