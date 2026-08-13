#!/usr/bin/env python3
"""
Build the candidate pool for negative DDIs from a Y2H/MS PPI parquet and dump it
to ``neg_pool.npz`` for the (seed-dependent) selection step.

This is the EXPENSIVE, DETERMINISTIC half of negative-DDI construction: it
streams the parquet, maps bait/prey genes to UniProt + Pfam via the UniProt REST
API, and assembles the pool of candidate Pfam pairs (restricted to Pfam domains
that already appear in a 3did positive DDI, and excluding pairs that already
exist as DDIs).  It performs NO sampling and NO insertion -- selection fans out
into parallel per-seed jobs (``select_ppi_negative_dans.py``) that read the dump,
and the winning selection is inserted by ``insert_ppi_negative_selection.py``.

The dump carries what both negative-construction methods need (uncapped DANS,
Cappelletti et al. vbae036):
  * Method 1 "deletion" -- the candidate pool ``cand_a``/``cand_b``, the pool
    domain universe ``pool_dom`` and the *reduced* positive degrees
    ``pool_deg_r`` (3did positives restricted to pool domains), plus the reduced
    positive-edge PA and target count.
  * Method 2 "random_addition" -- the full positive edge endpoint multiset
    ``pos_a``/``pos_b`` (DANS samples node-pairs proportional to degree from it),
    the full positive degrees/PA/count, and the forbidden-pair set
    ``forbidden_a``/``forbidden_b`` (all existing DDIs) that DANS must avoid.
"""

import argparse
import itertools
import json
import sqlite3
import sys
import time
import math
from collections import defaultdict

import numpy as np
import pyarrow.parquet as pq
import requests

from ddi_db_utils import pfam_sort_key


TAG = "[ppi_neg]"
BATCH_SIZE = 500_000
REQUIRED_COLUMNS = ["gene_name_bait", "gene_name_prey", "n_tested"]


def log(msg):
    print(f"{TAG} {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True)
    p.add_argument("--parquet", required=True)
    p.add_argument("--pfam-mapping-out", required=True,
                   help="Output path for UniProt -> Pfam JSON mapping")
    p.add_argument("--pool-out", required=True,
                   help="Output path for the candidate-pool .npz dump")
    p.add_argument("--min-n-tested", type=int, required=True)
    p.add_argument(
        "--no-self",
        action="store_true",
        help="Skip self-pairs (domain interacting with itself) "
             "if remove_self_interactions is enabled.",
    )
    return p.parse_args()


def fetch_gene_mappings(gene_names, batch_size=100):
    BASE_URL = "https://rest.uniprot.org/uniprotkb/search"
    HEADERS = {"accept": "application/json"}
    MAX_RETRIES = 10

    gene_list = sorted(gene_names)
    gene_to_uniprot = {}
    gene_seen = {}
    uniprot_to_pfams = defaultdict(set)
    n_batches = math.ceil(len(gene_list) / batch_size)

    log(f"fetching gene->UniProt+Pfam for {len(gene_list)} genes in {n_batches} batches")

    for i in range(0, len(gene_list), batch_size):
        batch = gene_list[i : i + batch_size]
        batch_num = i // batch_size + 1
        gene_clause = " OR ".join(f"gene:{g}" for g in batch)
        query = f"organism_name:Human AND ({gene_clause}) AND reviewed:true"
        params = {
            "query": query,
            "fields": "accession,gene_primary,xref_pfam",
            "sort": "accession desc",
            "size": "500",
        }

        for attempt in range(MAX_RETRIES):
            try:
                resp = requests.get(BASE_URL, headers=HEADERS, params=params, timeout=120)
                resp.raise_for_status()
                break
            except requests.RequestException as exc:
                if attempt < MAX_RETRIES - 1:
                    wait = min(2 ** (attempt + 1), 60)
                    log(f"batch {batch_num}/{n_batches} attempt {attempt + 1} failed: {exc}; retrying in {wait}s")
                    time.sleep(wait)
                else:
                    sys.exit(f"{TAG} batch {batch_num}/{n_batches} failed after {MAX_RETRIES} attempts: {exc}")

        data = resp.json()
        for entry in data.get("results", []):
            acc = entry.get("primaryAccession")
            genes = entry.get("genes", [])
            if not genes:
                continue
            primary_gene = genes[0].get("geneName", {}).get("value")
            if not primary_gene or primary_gene not in gene_names:
                for synonym in genes[0].get("synonyms", []):
                    syn_value = synonym.get("value")
                    if syn_value in gene_names:
                        primary_gene = syn_value
                        break
            if not primary_gene or primary_gene not in gene_names:
                continue

            if primary_gene not in gene_seen:
                gene_seen[primary_gene] = acc
            elif gene_seen[primary_gene] is not None and gene_seen[primary_gene] != acc:
                gene_seen[primary_gene] = None

            pfams = set()
            for xref in entry.get("uniProtKBCrossReferences", []):
                if xref.get("database") == "Pfam":
                    pfams.add(xref["id"].split(".")[0])
            if pfams:
                uniprot_to_pfams[acc] |= pfams

        log(f"batch {batch_num}/{n_batches}: queried {len(batch)} genes, got {len(data.get('results', []))} results")

        if batch_num < n_batches:
            time.sleep(1)

    n_ambig = 0
    ambiguous_genes = set()
    for gene, acc in gene_seen.items():
        if acc is not None:
            gene_to_uniprot[gene] = acc
        else:
            n_ambig += 1
            ambiguous_genes.add(gene)

    genes_not_seen = set(gene_names) - set(gene_seen.keys())

    log(f"n_mapped_genes = {len(gene_to_uniprot)}")
    log(f"n_unique_uniprots = {len(set(gene_to_uniprot.values()))}")
    log(f"n_dropped_ambiguous = {n_ambig}")
    log(",".join(ambiguous_genes))
    log(f"n_unseen_genes = {len(genes_not_seen)}")
    log(",".join(genes_not_seen))
    log(f"fetched Pfam mappings for {len(uniprot_to_pfams)} UniProt IDs")

    return gene_to_uniprot, uniprot_to_pfams


def load_3did_pfams(conn):
    """Pfam IDs that appear in a 3did positive DDI.

    Negatives are inferred (and degree-matched) only over the 3did domain
    universe, so single-domain / PPIDM positives never widen the candidate set.
    """
    cur = conn.execute(
        "SELECT DISTINCT d.pfam_id "
        "FROM domain AS d JOIN domain_domain_interaction AS ddi "
        "  ON d.id IN (ddi.domain_id_a, ddi.domain_id_b) "
        "WHERE ddi.negative = 0 AND ddi.source = '3did'"
    )
    return {row[0] for row in cur}


def load_existing_pairs(conn):
    cur = conn.execute(
        "SELECT da.pfam_id, db.pfam_id "
        "FROM domain_domain_interaction AS ddi "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b"
    )
    return {tuple(sorted((a, b), key=pfam_sort_key)) for a, b in cur}


def load_positive_3did_edges(conn):
    """The 3did positive DDIs as (pfam_a, pfam_b) pairs."""
    return conn.execute(
        "SELECT da.pfam_id, db.pfam_id "
        "FROM domain_domain_interaction AS ddi "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b "
        "WHERE ddi.negative = 0 AND ddi.source = '3did'"
    ).fetchall()


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

    gene_to_uniprot, uniprot_to_pfams = fetch_gene_mappings(unique_genes, batch_size=50)

    log(f"writing UniProt -> Pfam mapping to {args.pfam_mapping_out}")
    with open(args.pfam_mapping_out, "w") as fh:
        json.dump(
            {k: sorted(v) for k, v in uniprot_to_pfams.items()},
            fh,
        )

    n_pfam_unique = len({p for s in uniprot_to_pfams.values() for p in s})
    log(f"n_pfam_domains_for_input_proteins = {n_pfam_unique}")

    conn = sqlite3.connect(args.db)

    pos_pfam = load_3did_pfams(conn)
    log(f"n_3did_pfams = {len(pos_pfam)}")

    existing_pairs = load_existing_pairs(conn)
    log(f"n_existing_ddis = {len(existing_pairs)}")

    def row_pfams(gene):
        uniprot = gene_to_uniprot.get(gene)
        if uniprot is None:
            return set()
        return set(uniprot_to_pfams.get(uniprot, ()))

    candidate_counts = defaultdict(int)
    n_rows_with_pairs = 0
    n_skipped_unmapped = 0
    for bait, prey in zip(baits, preys):
        if bait not in gene_to_uniprot or prey not in gene_to_uniprot:
            n_skipped_unmapped += 1
            continue
        bait_pfams = row_pfams(bait) & pos_pfam
        prey_pfams = row_pfams(prey) & pos_pfam
        if not bait_pfams or not prey_pfams:
            continue
        row_pairs = set()
        for a, b in itertools.product(bait_pfams, prey_pfams):
            if args.no_self and a == b:
                continue
            row_pairs.add(tuple(sorted((a, b), key=pfam_sort_key)))
        if not row_pairs:
            continue
        n_rows_with_pairs += 1
        for key in row_pairs:
            candidate_counts[key] += 1

    del baits, preys

    log(f"n_ppi_rows_skipped_unmapped_gene = {n_skipped_unmapped}")
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

    # Positive 3did graph statistics: degree (the per-domain cap), edge PA, and
    # the target negative count.
    pos_edges = load_positive_3did_edges(conn)
    pos_degree = defaultdict(int)
    for a, b in pos_edges:
        pos_degree[a] += 1
        pos_degree[b] += 1
    n_positive = len(pos_edges)
    n_positive_domains = len(pos_degree)
    pos_edge_pa = np.array(
        [pos_degree[a] * pos_degree[b] for a, b in pos_edges], dtype=np.int64
    )
    log(f"n_positive_ddis_in_db = {n_positive}")
    log(f"n_positive_domains = {n_positive_domains}")
    log(f"positive mean PA = {float(pos_edge_pa.mean()):.1f}")

    # Drop candidates that already exist as a DDI (positive or negative).
    fresh_pairs = []
    n_positive_ddis_in_negative_ppis = 0
    for key in candidate_counts:
        if key in existing_pairs:
            n_positive_ddis_in_negative_ppis += 1
        else:
            fresh_pairs.append(key)

    log(f"n_positive_ddis_in_negative_ppis = {n_positive_ddis_in_negative_ppis}")
    log(f"n_fresh_candidates_after_dedup = {len(fresh_pairs)}")
    conn.close()

    cand_a = np.array([a for a, b in fresh_pairs], dtype=object)
    cand_b = np.array([b for a, b in fresh_pairs], dtype=object)

    # ---- Method 1 ("deletion"): reduce the positives to the candidate-domain
    #      universe so positive and candidate domains coincide; DANS then draws
    #      degree-aware over the fixed candidate pool. ----
    pool_domains = {d for pair in fresh_pairs for d in pair}
    pos_edges_r = [
        (a, b) for a, b in pos_edges if a in pool_domains and b in pool_domains
    ]
    pos_degree_r = defaultdict(int)
    for a, b in pos_edges_r:
        pos_degree_r[a] += 1
        pos_degree_r[b] += 1
    n_positive_r = len(pos_edges_r)
    n_positive_domains_r = len(pos_degree_r)
    pos_edge_pa_r = np.array(
        [pos_degree_r[a] * pos_degree_r[b] for a, b in pos_edges_r], dtype=np.int64
    )
    # Every pool domain carries its reduced-positive degree (0 if it has no edge
    # in the reduced positive graph); the selector turns these into PA weights.
    pool_dom = np.array(sorted(pool_domains, key=pfam_sort_key), dtype=object)
    pool_deg_r = np.array([pos_degree_r[d] for d in pool_dom], dtype=np.int64)
    log(f"n_pool_domains = {len(pool_dom)}")
    log(f"n_reduced_positive_ddis = {n_positive_r}")
    log(f"n_reduced_positive_domains = {n_positive_domains_r}")

    # ---- Method 2 ("random_addition"): plain DANS over the full positive set.
    #      The selector samples node-pairs proportional to degree by drawing from
    #      the endpoint multiset of these edges and rejects existing pairs. ----
    pos_a = np.array([a for a, b in pos_edges], dtype=object)
    pos_b = np.array([b for a, b in pos_edges], dtype=object)
    pos_dom = np.array(list(pos_degree.keys()), dtype=object)
    pos_deg = np.array([pos_degree[d] for d in pos_dom], dtype=np.int64)
    forbidden_a = np.array([a for a, b in existing_pairs], dtype=object)
    forbidden_b = np.array([b for a, b in existing_pairs], dtype=object)

    log(f"writing candidate pool to {args.pool_out}")
    np.savez(
        args.pool_out,
        # --- Method 1 (deletion) ---
        cand_a=cand_a,
        cand_b=cand_b,
        pool_dom=pool_dom,
        pool_deg_r=pool_deg_r,
        pos_edge_pa_r=pos_edge_pa_r,
        n_positive_r=np.int64(n_positive_r),
        n_positive_domains_r=np.int64(n_positive_domains_r),
        # --- Method 2 (random_addition) ---
        pos_a=pos_a,
        pos_b=pos_b,
        pos_dom=pos_dom,
        pos_deg=pos_deg,
        pos_edge_pa=pos_edge_pa,
        n_positive=np.int64(n_positive),
        n_positive_domains=np.int64(n_positive_domains),
        forbidden_a=forbidden_a,
        forbidden_b=forbidden_b,
    )
    log("done")


if __name__ == "__main__":
    main()
