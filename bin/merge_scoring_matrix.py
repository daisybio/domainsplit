#!/usr/bin/env python3
"""
merge_scoring_matrix.py
------------------------
Sums the raw, per-shard outputs of build_scoring_matrix.py (run with
--n_shards > 1) into the final, normalized scoring matrix for one method:

    <method>.c_ab_matrix.csv — 20x20 contact counts (aaA, aaB, count)
    <method>.db_freq.csv     — train-set surface AA frequencies (aa, frequency)
    <method>.t_db.txt        — total contact count T_DB

Output format is identical to the original (pre-sharding)
build_scoring_matrix.py, so score_ddi.py needs no changes to consume it.

Run once per method, after all of that method's build_scoring_matrix.py
shards have completed -- including the trivial --n_shards=1 case, so the
pipeline shape (shard -> merge) is uniform regardless of shard count.
"""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

from utils_struct import AA_3


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
    return p.parse_args()


def sum_c_ab_counts(paths):
    C_ab = defaultdict(lambda: defaultdict(int))
    for path in paths:
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                C_ab[row["aaA"]][row["aaB"]] += int(row["count"])
    return C_ab


def sum_surface_counts(paths):
    surface_counts = defaultdict(int)
    for path in paths:
        with open(path, newline="") as f:
            for row in csv.DictReader(f):
                surface_counts[row["aa"]] += int(row["count"])
    return surface_counts


def write_c_ab_matrix(C_ab, path):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["aaA", "aaB", "count"])
        for aaA in AA_3:
            for aaB in AA_3:
                writer.writerow([aaA, aaB, C_ab[aaA][aaB]])
    print(f"Written merged C_ab matrix to {path}", flush=True)


def write_db_freq(surface_counts, path):
    total = sum(surface_counts.values())
    if total == 0:
        print("WARNING: no surface residues accumulated — DB_FREQ will be zeros", file=sys.stderr)
        total = 1
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["aa", "frequency"])
        for aa in AA_3:
            writer.writerow([aa, surface_counts.get(aa, 0) / total])
    print(f"Written merged DB_FREQ to {path}", flush=True)


def write_t_db(C_ab, path):
    t_db_diag = int(sum(C_ab[aaA][aaA] for aaA in AA_3))
    t_db_offdiag = int(sum(
        C_ab[aaA][aaB] for i, aaA in enumerate(AA_3) for j, aaB in enumerate(AA_3) if j > i
    ))
    t_db = t_db_diag + 2 * t_db_offdiag
    Path(path).write_text(str(t_db) + "\n")
    print(f"T_DB = {t_db}, written to {path}", flush=True)
    return t_db


def main():
    args = parse_args()
    print(f"[merge_scoring_matrix] method={args.method}: merging "
          f"{len(args.c_ab_counts)} shard(s)", flush=True)

    C_ab = sum_c_ab_counts(args.c_ab_counts)
    surface_counts = sum_surface_counts(args.surface_counts)

    write_c_ab_matrix(C_ab, args.c_ab_matrix_out)
    write_db_freq(surface_counts, args.db_freq_out)
    t_db = write_t_db(C_ab, args.t_db_out)

    print(f"[merge_scoring_matrix] method={args.method}: done, T_DB={t_db}", flush=True)


if __name__ == "__main__":
    main()
