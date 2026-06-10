#!/usr/bin/env python3
"""
Insert the two negative-DDI construction methods into the domainsplit SQLite,
each under its own ``source`` labels so both can coexist in one database and be
selected independently by the downstream splits.

Both methods use uncapped Degree-Aware Node Sampling (DANS); the per-seed pairs
are produced by ``select_ppi_negative_dans.py`` (one run per method, no
pick-best).  This step copies the matching positives and inserts the negatives:

  * "deletion"        -- positives = 3did restricted to the candidate-pool domain
    universe (``3did_deletion``); negatives =
    ``inferred_ppi_screen_negative_for_deletion``.
  * "random_addition" -- positives = the full 3did set (``3did_random_addition``);
    negatives = ``inferred_ppi_screen_negative_for_random_addition``.

The four method labels are inserted with ``dedup_across_sources=False`` so they
may duplicate a pair already stored under the canonical ``3did`` source (the
table's ``UNIQUE(domain_id_a, domain_id_b, source)`` keeps the labels distinct).

Prints a positive reference + the negative selection for each method and writes
the same data to a published scores TSV.
"""

import argparse
import json
import sqlite3

import numpy as np

from ddi_db_utils import count_source, ensure_domains, insert_ddis


TAG = "[neg_insert]"


def log(msg):
    print(f"{TAG} {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True)
    p.add_argument("--pool", required=True,
                   help="candidate-pool .npz (for the pool-domain universe)")
    p.add_argument("--pairs-deletion", required=True)
    p.add_argument("--pairs-random-addition", required=True)
    p.add_argument("--score-deletion", required=True)
    p.add_argument("--score-random-addition", required=True)
    p.add_argument("--scores-out", required=True,
                   help="output consolidated scores TSV path")
    p.add_argument("--source-3did", default="3did")
    p.add_argument("--label-pos-deletion", default="3did_deletion")
    p.add_argument("--label-pos-random-addition", default="3did_random_addition")
    p.add_argument("--label-neg-deletion",
                   default="inferred_ppi_screen_negative_for_deletion")
    p.add_argument("--label-neg-random-addition",
                   default="inferred_ppi_screen_negative_for_random_addition")
    return p.parse_args()


def read_pairs(path):
    pairs = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            a, b = line.split("\t")
            pairs.append((a, b))
    return pairs


def load_positive_pairs(conn, source):
    """3did positive (pfam_a, pfam_b) pairs currently in the DB."""
    return conn.execute(
        "SELECT da.pfam_id, db.pfam_id "
        "FROM domain_domain_interaction AS ddi "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b "
        "WHERE ddi.negative = 0 AND ddi.source = ?",
        (source,),
    ).fetchall()


def main():
    args = parse_args()

    pool_domains = set(np.load(args.pool, allow_pickle=True)["pool_dom"].tolist())
    log(f"pool-domain universe: {len(pool_domains)} domains")

    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    # --- positives: full (random_addition) and pool-restricted (deletion) copies ---
    pos_3did = load_positive_pairs(conn, args.source_3did)
    pos_reduced = [
        (a, b) for a, b in pos_3did if a in pool_domains and b in pool_domains
    ]
    log(f"3did positives: {len(pos_3did)} total, {len(pos_reduced)} within pool")

    insert_ddis(conn, pos_3did, negative=False,
                source=args.label_pos_random_addition, dedup_across_sources=False)
    insert_ddis(conn, pos_reduced, negative=False,
                source=args.label_pos_deletion, dedup_across_sources=False)

    # --- negatives: one DANS selection per method ---
    pairs_del = read_pairs(args.pairs_deletion)
    pairs_rand = read_pairs(args.pairs_random_addition)

    ensure_domains(conn, (p for pair in pairs_del for p in pair))
    ensure_domains(conn, (p for pair in pairs_rand for p in pair))
    insert_ddis(conn, pairs_del, negative=True,
                source=args.label_neg_deletion, dedup_across_sources=False)
    insert_ddis(conn, pairs_rand, negative=True,
                source=args.label_neg_random_addition, dedup_across_sources=False)

    conn.commit()
    counts = {
        args.label_pos_random_addition: count_source(conn, args.label_pos_random_addition),
        args.label_pos_deletion: count_source(conn, args.label_pos_deletion),
        args.label_neg_random_addition: count_source(conn, args.label_neg_random_addition),
        args.label_neg_deletion: count_source(conn, args.label_neg_deletion),
    }
    conn.close()

    with open(args.score_deletion) as fh:
        sc_del = json.load(fh)
    with open(args.score_random_addition) as fh:
        sc_rand = json.load(fh)

    for label, n in counts.items():
        log(f"inserted source={label} rows={n}")
    for sc in (sc_del, sc_rand):
        log(f"method={sc['method']} set=positive n_sel={sc['pos_n_sel']} "
            f"n_dom={sc['pos_n_dom']} mean_pa={sc['pos_mean_pa']:.1f}")
        log(f"method={sc['method']} set=negative seed={sc['seed']} J={sc['J']:.4f} "
            f"pa={sc['pa']:.4f} deg={sc['deg']:.4f} cov={sc['cov']:.4f} "
            f"n_sel={sc['n_sel']} n_dom={sc['n_dom']} mean_pa={sc['mean_pa']:.1f}")

    cols = ["set", "method", "seed", "J", "pa", "deg", "cov", "n_sel", "n_dom",
            "mean_pa"]
    with open(args.scores_out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for sc in (sc_del, sc_rand):
            fh.write("\t".join([
                "positive", sc["method"], "NA", "NA", "NA", "NA", "NA",
                str(sc["pos_n_sel"]), str(sc["pos_n_dom"]),
                f"{sc['pos_mean_pa']:.2f}",
            ]) + "\n")
            fh.write("\t".join([
                "negative", sc["method"], str(sc["seed"]),
                f"{sc['J']:.6f}", f"{sc['pa']:.6f}", f"{sc['deg']:.6f}",
                f"{sc['cov']:.6f}", str(sc["n_sel"]), str(sc["n_dom"]),
                f"{sc['mean_pa']:.2f}",
            ]) + "\n")


if __name__ == "__main__":
    main()
