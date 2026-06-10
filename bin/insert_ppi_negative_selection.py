#!/usr/bin/env python3
"""
Pick the best per-seed negative-DDI selection and insert it into the domainsplit
SQLite.

Reads every ``score_*.json`` + ``pairs_*.tsv`` produced by the parallel
``select_ppi_negative_dans.py`` jobs, picks the selection with the lowest
objective ``J`` (ties broken by the smaller seed, so the result is fully
deterministic regardless of which SLURM task finished first), and inserts that
seed's Pfam pairs as negatives via the shared ``ddi_db_utils`` helpers.

Prints a positive reference line (absolute baseline) followed by one line per
seed and a WINNER line, and writes the same data to a published scores TSV.
"""

import argparse
import glob
import json
import sqlite3

from ddi_db_utils import count_source, ensure_domains, insert_ddis


TAG = "[neg_insert]"
J_ROUND = 12


def log(msg):
    print(f"{TAG} {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True)
    p.add_argument("--scores-out", required=True,
                   help="output consolidated scores TSV path")
    p.add_argument("--score-glob", default="score_*.json")
    p.add_argument("--pairs-template", default="pairs_{seed}.tsv")
    p.add_argument("--source-label", default="inferred_ppi_screen_negative")
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


def main():
    args = parse_args()

    score_files = sorted(glob.glob(args.score_glob))
    if not score_files:
        raise SystemExit(f"{TAG} no score files matching {args.score_glob}")

    records = []
    for path in score_files:
        with open(path) as fh:
            records.append(json.load(fh))

    # Deterministic winner: lowest J (rounded), ties broken by smaller seed.
    records.sort(key=lambda r: (round(r["J"], J_ROUND), r["seed"]))
    winner = records[0]
    winner_seed = winner["seed"]

    pairs_path = args.pairs_template.format(seed=winner_seed)
    pairs = read_pairs(pairs_path)

    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    ensure_domains(conn, (p for pair in pairs for p in pair))
    insert_ddis(conn, pairs, negative=True, source=args.source_label)
    conn.commit()
    n_inserted = count_source(conn, args.source_label)
    conn.close()

    # --- report: positive reference, then every seed, then the winner ---
    ref = records[0]
    log(f"set=positive n_sel={ref['pos_n_sel']} n_dom={ref['pos_n_dom']} "
        f"mean_pa={ref['pos_mean_pa']:.1f}")
    for r in records:
        log(f"set=negative seed={r['seed']} J={r['J']:.4f} pa={r['pa']:.4f} "
            f"deg={r['deg']:.4f} cov={r['cov']:.4f} n_sel={r['n_sel']} "
            f"n_dom={r['n_dom']} mean_pa={r['mean_pa']:.1f}")
    log(f"WINNER seed={winner_seed} J={winner['J']:.4f} n_inserted={n_inserted}")

    cols = ["set", "seed", "J", "pa", "deg", "cov", "n_sel", "n_dom",
            "mean_pa", "winner"]
    with open(args.scores_out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        fh.write("\t".join([
            "positive", "NA", "NA", "NA", "NA", "NA",
            str(ref["pos_n_sel"]), str(ref["pos_n_dom"]),
            f"{ref['pos_mean_pa']:.2f}", "NA",
        ]) + "\n")
        for r in records:
            fh.write("\t".join([
                "negative", str(r["seed"]),
                f"{r['J']:.6f}", f"{r['pa']:.6f}", f"{r['deg']:.6f}",
                f"{r['cov']:.6f}", str(r["n_sel"]), str(r["n_dom"]),
                f"{r['mean_pa']:.2f}",
                "1" if r["seed"] == winner_seed else "0",
            ]) + "\n")


if __name__ == "__main__":
    main()
