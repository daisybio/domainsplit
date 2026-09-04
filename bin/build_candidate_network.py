#!/usr/bin/env python3
"""Build the candidate network of high-confidence non-interacting domain pairs.

Streams a Y2H/MS PPI parquet, maps bait/prey genes to UniProt + Pfam via the
UniProt REST API, and writes every Pfam pair the screen tested without finding an
interaction as ``candidate_network.csv``.  ppi-splitting consumes that file as
the pool its ``ilp_candidates`` negative sampler draws from -- these are
*experimentally supported* non-interactions, as opposed to the uniform sampler's
"any pair not known to interact".

No sampling happens here: the file is the full pool, and ppi-splitting decides
how much of it to use per split.

Two restrictions define the pool:
  * both families must appear in a **3did** DDI -- 3did is the population being
    partitioned, so a negative over any other family universe could never be
    paired against a positive;
  * pairs that already exist as a DDI are excluded.  At the point this runs the
    database holds 3did only (the external sources are still TSVs, inserted
    after ppi-splitting returns), so this excludes exactly the 3did positives.

Self-pairs (a family against itself) are kept: they are legitimate DDIs, and
ppi-splitting represents them as instance pairs drawn from the same family.

The gene -> UniProt -> Pfam lookup is the only network call in this pipeline that
cannot be replaced by a downloaded file, so it is cacheable: the mapping this
script writes (``--mapping-out``) can be handed straight back to it
(``--mapping-in``), and then no request is made at all. That is what makes the
`-profile test` run offline, and it also stops a re-run of a real analysis from
re-querying tens of thousands of genes.
"""

import argparse
import csv
import itertools
import json
import math
import sqlite3
import sys
import time
from collections import defaultdict

import pyarrow.parquet as pq
import requests

from ddi_db_utils import canonical_pair, pfam_sort_key

TAG = "[candidate_network]"
BATCH_SIZE = 500_000
REQUIRED_COLUMNS = ["gene_name_bait", "gene_name_prey", "n_tested"]


def log(msg):
    print(f"{TAG} {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True)
    p.add_argument("--parquet", required=True)
    p.add_argument("--mapping-out", required=True,
                   help="Output path for the gene -> UniProt -> Pfam JSON mapping")
    p.add_argument("--mapping-in",
                   help="A previous --mapping-out file. Given, no UniProt request is made; "
                        "genes it does not cover are treated as unmapped.")
    p.add_argument("--network-out", required=True,
                   help="Output path for the candidate_network CSV (protein1,protein2)")
    p.add_argument("--min-n-tested", type=int, required=True)
    return p.parse_args()


def read_mapping(path):
    """``(gene_to_uniprot, uniprot_to_pfams)`` from a previous ``--mapping-out``."""
    with open(path) as fh:
        data = json.load(fh)
    missing = {"gene_to_uniprot", "uniprot_to_pfams"} - set(data)
    if missing:
        raise ValueError(
            f"{path} is not a candidate-network mapping file (missing {sorted(missing)})"
        )
    return (
        dict(data["gene_to_uniprot"]),
        {k: set(v) for k, v in data["uniprot_to_pfams"].items()},
    )


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
    log(",".join(sorted(ambiguous_genes)))
    log(f"n_unseen_genes = {len(genes_not_seen)}")
    log(",".join(sorted(genes_not_seen)))
    log(f"fetched Pfam mappings for {len(uniprot_to_pfams)} UniProt IDs")

    return gene_to_uniprot, uniprot_to_pfams


def load_3did_pfams(conn):
    """Pfam IDs that appear in a 3did positive DDI.

    Candidates are restricted to this universe so a sampled negative always sits
    in the same family population as the positives it is paired against.
    """
    cur = conn.execute(
        "SELECT DISTINCT d.pfam_id "
        "FROM domain AS d JOIN domain_domain_interaction AS ddi "
        "  ON d.id IN (ddi.domain_id_a, ddi.domain_id_b) "
        "WHERE ddi.negative = 0 "
        "  AND ',' || ddi.source || ',' LIKE '%,3did,%'"
    )
    return {row[0] for row in cur}


def load_existing_pairs(conn):
    cur = conn.execute(
        "SELECT da.pfam_id, db.pfam_id "
        "FROM domain_domain_interaction AS ddi "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b"
    )
    return {canonical_pair(a, b) for a, b in cur}


def _validate_columns(parquet_schema):
    available = set(parquet_schema.names)
    missing = set(REQUIRED_COLUMNS) - available
    if missing:
        sys.exit(f"{TAG} parquet missing required columns: {missing}")


def _collect_genes_and_pairs(parquet_path, min_n_tested):
    """Stream parquet in batches. Returns (n_input, unique_genes, baits, preys)."""
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
    log(f"n_input_ppis = {n_input}")
    log(f"n_ppis_after_n_tested_filter (>= {args.min_n_tested}) = {len(baits)}")
    log(f"n_unique_genes = {len(unique_genes)}")

    if args.mapping_in:
        gene_to_uniprot, uniprot_to_pfams = read_mapping(args.mapping_in)
        covered = len(unique_genes & set(gene_to_uniprot))
        log(f"reusing mapping from {args.mapping_in}: no UniProt request; "
            f"{covered}/{len(unique_genes)} parquet genes covered")
    else:
        gene_to_uniprot, uniprot_to_pfams = fetch_gene_mappings(unique_genes, batch_size=50)

    log(f"writing gene -> UniProt -> Pfam mapping to {args.mapping_out}")
    with open(args.mapping_out, "w") as fh:
        json.dump(
            {
                "gene_to_uniprot": dict(sorted(gene_to_uniprot.items())),
                "uniprot_to_pfams": {k: sorted(v) for k, v in sorted(uniprot_to_pfams.items())},
            },
            fh,
            indent=1,
        )

    n_pfam_unique = len({p for s in uniprot_to_pfams.values() for p in s})
    log(f"n_pfam_domains_for_input_proteins = {n_pfam_unique}")

    conn = sqlite3.connect(args.db)
    pos_pfam = load_3did_pfams(conn)
    log(f"n_3did_pfams = {len(pos_pfam)}")
    existing_pairs = load_existing_pairs(conn)
    log(f"n_existing_ddis = {len(existing_pairs)}")
    conn.close()

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
        row_pairs = {
            canonical_pair(a, b)
            for a, b in itertools.product(bait_pfams, prey_pfams)
        }
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
        log(f"most_common_ddi = {most_common_key[0]}-{most_common_key[1]} "
            f"(observed in {most_common_count} PPI rows)")

    fresh_pairs = [k for k in candidate_counts if k not in existing_pairs]
    log(f"n_candidates_dropped_as_existing_ddi = "
        f"{len(candidate_counts) - len(fresh_pairs)}")
    log(f"n_candidate_network_pairs = {len(fresh_pairs)}")

    fresh_pairs.sort(key=lambda pair: (pfam_sort_key(pair[0]), pfam_sort_key(pair[1])))
    with open(args.network_out, "w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["protein1", "protein2"])
        writer.writerows(fresh_pairs)
    log(f"wrote {args.network_out}")


if __name__ == "__main__":
    main()
