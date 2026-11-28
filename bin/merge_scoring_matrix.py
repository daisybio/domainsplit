#!/usr/bin/env python3
"""
merge_scoring_matrix.py
------------------------
Sums the raw, per-shard outputs of build_scoring_matrix.py (run with
--n_shards > 1) into the final, normalized scoring matrix for one method:

    <method>.c_ab_matrix.csv — 20x20 contact counts (aaA, aaB, count)
    <method>.db_freq.csv     — train-set surface AA frequencies (aa, frequency)
    <method>.t_db.txt        — total contact count T_DB
    <method>.fold_counts.json — {"c_ab": [fold][aaA][aaB], "surface": [fold][aa]}
                                raw counts per cross-fit fold, AA_3 order;
                                score_ddi.py subtracts a train row's own fold

The first three are the sum over all folds.

Run once per method, after all of that method's build_scoring_matrix.py
shards have completed -- including the trivial --n_shards=1 case, so the
pipeline shape (shard -> merge) is uniform regardless of shard count.
"""
import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

from utils_struct import AA_3, N_FOLDS


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--method", required=True)
    p.add_argument("--c_ab_counts", required=True, nargs="+",
                    help="one or more <method>.<shard>.c_ab_counts.csv files to sum")
    p.add_argument("--surface_counts", required=True, nargs="+",
                    help="one or more <method>.<shard>.surface_counts.csv files to sum")
    p.add_argument("--c_ab_matrix_out", required=True)
    p.add_argument("--db_freq_out", required=True)
    p.add_argument("--t_db_out", required=True)
    p.add_argument("--fold_counts_out", required=True)
    return p.parse_args()


IDX = {aa: i for i, aa in enumerate(AA_3)}


def sum_c_ab_counts(paths):
    """(fold, aaA, aaB) counts summed over shards, AA_3 order."""
    C_ab = np.zeros((N_FOLDS, len(AA_3), len(AA_3)), dtype=np.int64)
    for path in paths:
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                C_ab[int(row["fold"]), IDX[row["aaA"]], IDX[row["aaB"]]] += int(row["count"])
    return C_ab


def sum_surface_counts(paths):
    """(fold, aa) counts summed over shards, AA_3 order."""
    surface_counts = np.zeros((N_FOLDS, len(AA_3)), dtype=np.int64)
    for path in paths:
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                surface_counts[int(row["fold"]), IDX[row["aa"]]] += int(row["count"])
    return surface_counts


def write_fold_counts(C_ab, surface_counts, path):
    with open(path, "w") as f:
        json.dump({"c_ab": C_ab.tolist(), "surface": surface_counts.tolist()}, f)
    print(f"Written per-fold counts to {path}", flush=True)


def write_c_ab_matrix(C_ab, path):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["aaA", "aaB", "count"])
        for i, aaA in enumerate(AA_3):
            for j, aaB in enumerate(AA_3):
                writer.writerow([aaA, aaB, int(C_ab[i, j])])
    print(f"Written merged C_ab matrix to {path}", flush=True)


def write_db_freq(surface_counts, path):
    total = int(surface_counts.sum())
    if total == 0:
        print("WARNING: no surface residues accumulated — DB_FREQ will be zeros", file=sys.stderr)
        total = 1
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["aa", "frequency"])
        for i, aa in enumerate(AA_3):
            writer.writerow([aa, int(surface_counts[i]) / total])
    print(f"Written merged DB_FREQ to {path}", flush=True)


def write_t_db(C_ab, path):
    t_db = int(np.trace(C_ab) + 2 * np.triu(C_ab, 1).sum())
    Path(path).write_text(str(t_db) + "\n")
    print(f"T_DB = {t_db}, written to {path}", flush=True)
    return t_db


def main():
    args = parse_args()
    print(f"[merge_scoring_matrix] method={args.method}: merging "
          f"{len(args.c_ab_counts)} shard(s)", flush=True)

    C_ab_folds = sum_c_ab_counts(args.c_ab_counts)
    surface_folds = sum_surface_counts(args.surface_counts)
    write_fold_counts(C_ab_folds, surface_folds, args.fold_counts_out)

    C_ab = C_ab_folds.sum(0)
    write_c_ab_matrix(C_ab, args.c_ab_matrix_out)
    write_db_freq(surface_folds.sum(0), args.db_freq_out)
    t_db = write_t_db(C_ab, args.t_db_out)

    print(f"[merge_scoring_matrix] method={args.method}: done, T_DB={t_db}", flush=True)


if __name__ == "__main__":
    main()
