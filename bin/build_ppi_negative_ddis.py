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
import re
import sqlite3
import sys
import time
from collections import defaultdict

import math

import numpy as np
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
    p.add_argument("--swissprot", required=True, help="UniProt Swiss-Prot .dat or .dat.gz")
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


def fetch_pfam_from_uniprot(uniprot_ids, batch_size=500):
    BASE_URL = "https://rest.uniprot.org/uniprotkb/search"
    HEADERS = {"accept": "application/json"}
    MAX_RETRIES = 3

    id_list = sorted(uniprot_ids)
    uniprot_to_pfams = defaultdict(set)
    n_batches = math.ceil(len(id_list) / batch_size)

    log(f"fetching Pfam cross-refs for {len(id_list)} UniProt IDs in {n_batches} batches")

    for i in range(0, len(id_list), batch_size):
        batch = id_list[i : i + batch_size]
        batch_num = i // batch_size + 1
        query = " OR ".join(batch)
        params = {
            "query": query,
            "fields": "accession,xref_pfam",
            "size": str(batch_size),
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
            if acc not in uniprot_ids:
                continue
            for xref in entry.get("uniProtKBCrossReferences", []):
                if xref.get("database") == "Pfam":
                    pfam_id = xref["id"].split(".")[0]
                    uniprot_to_pfams[acc].add(pfam_id)

        log(f"batch {batch_num}/{n_batches}: queried {len(batch)} IDs, got {len(data.get('results', []))} results")

    log(f"fetched Pfam mappings for {len(uniprot_to_pfams)} UniProt IDs")
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


def _strip_evidence(s):
    return s.split("{")[0].strip()


def parse_swissprot(path, gene_set):
    """Parse uniprot_sprot.dat(.gz) into gene->uniprot mapping (primary AC only).

    Name-level mappings take priority over Synonym-level.
    Genes mapping to >1 UniProt entry at the same level are dropped.
    """
    gene_name_map = {}
    gene_synonym_map = {}
    n_entries = 0

    with open_idmapping(path) as fh:
        primary_ac = None
        gn_lines = []

        for line in fh:
            if line.startswith("AC   "):
                if primary_ac is None:
                    primary_ac = line[5:].strip().split(";")[0].strip()
            elif line.startswith("GN   "):
                gn_lines.append(line[5:].rstrip("\n"))
            elif line.startswith("//"):
                if primary_ac and gn_lines:
                    n_entries += 1
                    gn_text = " ".join(gn_lines)

                    names = [_strip_evidence(n) for n in re.findall(r"Name=([^;]+);", gn_text)]
                    synonyms = []
                    for syn_group in re.findall(r"Synonyms=([^;]+);", gn_text):
                        synonyms.extend(_strip_evidence(s) for s in syn_group.split(","))

                    for gene in names:
                        if gene not in gene_set:
                            continue
                        if gene not in gene_name_map:
                            gene_name_map[gene] = primary_ac
                        elif gene_name_map[gene] is not None and gene_name_map[gene] != primary_ac:
                            gene_name_map[gene] = None

                    for gene in synonyms:
                        if gene not in gene_set:
                            continue
                        if gene not in gene_synonym_map:
                            gene_synonym_map[gene] = primary_ac
                        elif gene_synonym_map[gene] is not None and gene_synonym_map[gene] != primary_ac:
                            gene_synonym_map[gene] = None

                primary_ac = None
                gn_lines = []

    gene_to_uniprot = {}
    n_ambig_name = 0
    n_ambig_synonym = 0

    for gene in set(gene_name_map) | set(gene_synonym_map):
        name_ac = gene_name_map.get(gene)
        if gene in gene_name_map:
            if name_ac is not None:
                gene_to_uniprot[gene] = name_ac
            else:
                n_ambig_name += 1
        else:
            syn_ac = gene_synonym_map.get(gene)
            if syn_ac is not None:
                gene_to_uniprot[gene] = syn_ac
            else:
                n_ambig_synonym += 1

    n_dropped = n_ambig_name + n_ambig_synonym
    log(f"parsed {n_entries} SwissProt entries with gene names")
    log(f"n_mapped_genes = {len(gene_to_uniprot)}")
    log(f"n_unique_uniprots = {len(set(gene_to_uniprot.values()))}")
    log(f"n_dropped_ambiguous = {n_dropped} (name={n_ambig_name}, synonym={n_ambig_synonym})")

    return gene_to_uniprot


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

    log(f"parsing SwissProt file {args.swissprot}")
    gene_to_uniprot = parse_swissprot(args.swissprot, unique_genes)

    query_uniprots = set(gene_to_uniprot.values())
    log(f"n_unique_uniprots_to_query = {len(query_uniprots)}")
    uniprot_to_pfams = fetch_pfam_from_uniprot(query_uniprots)

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
            if a == b:
                continue
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
