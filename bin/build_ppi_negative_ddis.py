#!/usr/bin/env python3
"""
Build negative DDIs from a Y2H/MS PPI parquet and append them to the
domainsplit SQLite, restricted to Pfam domains already present in positive
DDIs.  In degree_matched mode, selects pairs so each domain's negative
degree matches its positive degree.  In frequency mode, takes the top-N
by PPI co-occurrence count, capped at (n_positive - n_negatome).
"""

import argparse
import heapq
import itertools
import json
import random
import sqlite3
import sys
import time
import math
from collections import defaultdict

import pyarrow.parquet as pq
import requests


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


def fetch_gene_mappings(gene_names, batch_size=100):
    BASE_URL = "https://rest.uniprot.org/uniprotkb/search"
    HEADERS = {"accept": "application/json"}
    MAX_RETRIES = 3

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
                    wait = 2 ** (attempt + 1)
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
    """Per-Pfam degree in the positive DDI set."""
    rows = conn.execute(
        "SELECT da.pfam_id, db.pfam_id "
        "FROM domain_domain_interaction AS ddi "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b "
        "WHERE ddi.negative = 0"
    ).fetchall()
    deg = defaultdict(int)
    for a, b in rows:
        deg[a] += 1
        deg[b] += 1
    return deg


def select_degree_matched(fresh_candidates, pos_degree, n_take):
    """Select negatives so each domain's negative degree matches its positive degree.

    Uses a lazy-deletion max-heap scored by combined degree deficit of both
    domains in each candidate pair.  Candidates are shuffled for random
    tiebreaking among equal-deficit pairs.
    """
    if not fresh_candidates or n_take <= 0:
        return []
    if not pos_degree:
        return fresh_candidates[:n_take]

    target = dict(pos_degree)
    current = defaultdict(int)

    candidates = list(fresh_candidates)
    random.shuffle(candidates)

    remaining = set(range(len(candidates)))

    def deficit(pfam):
        return max(0, target.get(pfam, 0) - current[pfam])

    def score(i):
        (pfam_a, pfam_b), _ = candidates[i]
        return deficit(pfam_a) + deficit(pfam_b)

    heap = [(-score(i), i) for i in range(len(candidates))]
    heapq.heapify(heap)

    chosen = []
    while len(chosen) < n_take and heap:
        neg_s, i = heapq.heappop(heap)
        if i not in remaining:
            continue

        actual = score(i)
        if actual != -neg_s:
            if actual > 0:
                heapq.heappush(heap, (-actual, i))
            else:
                remaining.discard(i)
            continue

        if actual <= 0:
            break

        (pfam_a, pfam_b), count = candidates[i]
        chosen.append(((pfam_a, pfam_b), count))
        remaining.discard(i)
        current[pfam_a] += 1
        current[pfam_b] += 1

    matched = sum(1 for p in target if current.get(p, 0) >= target[p])
    over = sum(1 for p in target if current.get(p, 0) > target[p])
    total_deficit = sum(max(0, target[p] - current.get(p, 0)) for p in target)
    log(f"degree_matched: {matched}/{len(target)} domains reached target degree")
    log(f"degree_matched: {over} domains exceeded target degree")
    log(f"degree_matched: remaining total deficit = {total_deficit}")
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

    gene_to_uniprot, uniprot_to_pfams = fetch_gene_mappings(unique_genes)

    log(f"writing UniProt -> Pfam mapping to {args.pfam_mapping_out}")
    with open(args.pfam_mapping_out, "w") as fh:
        json.dump(
            {k: sorted(v) for k, v in uniprot_to_pfams.items()},
            fh,
        )

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
            row_pairs.add(tuple(sorted((a, b))))
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

    n_positive = conn.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction WHERE negative = 0"
    ).fetchone()[0]
    n_negatome = conn.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction "
        "WHERE negative = 1 AND source = 'negatome'"
    ).fetchone()[0]
    log(f"n_positive_ddis_in_db = {n_positive}")
    log(f"n_negatome_negatives_in_db = {n_negatome}")
    if args.sampling_strategy == "degree_matched":
        n_take = n_positive
        log(f"n_take = {n_take} (degree_matched: matching positive count)")
    else:
        n_take = max(0, n_positive - n_negatome)
        log(f"n_take (target for source='{args.source_label}') = {n_take}")

    fresh_candidates = []
    n_positive_ddis_in_negative_ppis = 0

    for key, count in candidate_counts.items():
        if key in existing_pairs:
            n_positive_ddis_in_negative_ppis += 1
        else:
            fresh_candidates.append((key, count))

    log(f"n_positive_ddis_in_negative_ppis = {n_positive_ddis_in_negative_ppis}")

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
        # Pre-load pfam_id -> domain.id mapping to avoid per-row subqueries
        pfam_to_domain_ids = defaultdict(list)
        for did, pfam in conn.execute("SELECT id, pfam_id FROM domain"):
            pfam_to_domain_ids[pfam].append(did)
        log(f"loaded {len(pfam_to_domain_ids)} pfam -> domain mappings")

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
