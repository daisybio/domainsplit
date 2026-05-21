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
import sqlite3
import sys
from collections import defaultdict

import pandas as pd


TAG = "[ppi_neg]"


def log(msg):
    print(f"{TAG} {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True)
    p.add_argument("--parquet", required=True)
    p.add_argument("--idmapping", required=True, help="UniProt idmapping .dat or .dat.gz")
    p.add_argument("--min-n-tested", type=int, required=True)
    p.add_argument("--source-label", required=True)
    return p.parse_args()


def open_idmapping(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


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
    """One pass over idmapping. Returns (gene_to_uniprots, uniprot_to_pfams).

    Only keeps Gene_Name / Gene_Synonym rows whose value is in ``gene_set``,
    and only keeps Pfam rows for the UniProts that survive that filter.
    """
    gene_to_uniprots = defaultdict(set)
    relevant_uniprots = set()
    pfam_buffer = []
    with open_idmapping(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            uniprot, kind, value = parts
            if kind == "Gene_Name" or kind == "Gene_Synonym":
                if value in gene_set:
                    gene_to_uniprots[value].add(uniprot)
                    relevant_uniprots.add(uniprot)
            elif kind == "Pfam":
                pfam_buffer.append((uniprot, value))

    uniprot_to_pfams = defaultdict(set)
    for uniprot, pfam in pfam_buffer:
        if uniprot in relevant_uniprots:
            uniprot_to_pfams[uniprot].add(pfam.split(".")[0])
    return gene_to_uniprots, uniprot_to_pfams


def main():
    args = parse_args()

    log(f"loading parquet {args.parquet}")
    df = pd.read_parquet(args.parquet)
    n_input = len(df)
    log(f"n_input_ppis = {n_input}")

    required = {"gene_name_bait", "gene_name_prey", "n_tested"}
    missing = required - set(df.columns)
    if missing:
        sys.exit(f"{TAG} parquet missing required columns: {missing}")

    df = df[df["n_tested"] >= args.min_n_tested]
    df = df.dropna(subset=["gene_name_bait", "gene_name_prey"])
    n_after = len(df)
    log(f"n_ppis_after_n_tested_filter (>= {args.min_n_tested}) = {n_after}")

    unique_genes = set(df["gene_name_bait"]).union(df["gene_name_prey"])
    log(f"n_unique_genes = {len(unique_genes)}")

    log(f"streaming idmapping {args.idmapping}")
    gene_to_uniprots, uniprot_to_pfams = stream_idmapping(args.idmapping, unique_genes)
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
    for bait, prey in zip(df["gene_name_bait"], df["gene_name_prey"]):
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
    chosen = fresh_candidates[:n_take]
    log(f"n_chosen = {len(chosen)}")

    if chosen:
        conn.executemany(
            "INSERT OR IGNORE INTO domain_domain_interaction"
            "(domain_id_a, domain_id_b, negative, source) "
            "SELECT d1.id, d2.id, TRUE, ? "
            "FROM domain AS d1, domain AS d2 "
            "WHERE d1.pfam_id = ? AND d2.pfam_id = ?;",
            [(args.source_label, a, b) for (a, b), _ in chosen],
        )
        conn.commit()

    n_inserted = conn.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction WHERE source = ?",
        (args.source_label,),
    ).fetchone()[0]
    log(f"n_inserted_for_source = {n_inserted}")
    conn.close()


if __name__ == "__main__":
    main()
