#!/usr/bin/env python3
"""Merge every method's scores.tsv/confirmed.tsv (score_ddi.py) into one copy
of the struct-enriched master. Each method's files touch only its own
ddi_split_membership rows and ddi_interaction_confirmed rows, so the merges
are disjoint and order-independent."""
import argparse as ap
import os
import shutil
import sqlite3
import sys


def parse_args():
    p = ap.ArgumentParser()
    p.add_argument("--db_in", required=True)
    p.add_argument("--db_out", required=True)
    p.add_argument("--scores", required=True, nargs="+")
    p.add_argument("--confirmed", required=True, nargs="+")
    p.add_argument("--versions", required=True)
    p.add_argument("--process_name", required=True)
    return p.parse_args()


def apply_scores(conn, tsv_path):
    n = 0
    with open(tsv_path) as fh:
        next(fh)
        for line in fh:
            ddi_id, method, ia, ib, z_score = line.rstrip("\n").split("\t")
            # Every split row of the triple: the 1:1 and 1:10 test splits share
            # their positives, scored once in score_ddi.py.
            conn.execute(
                "UPDATE ddi_split_membership SET z_score = ? "
                "WHERE ddi_id = ? AND method = ? "
                "AND instance_id_a = ? AND instance_id_b = ?",
                (float(z_score), int(ddi_id), method, ia, ib),
            )
            n += 1
    return n


def apply_confirmed(conn, tsv_path):
    n = 0
    with open(tsv_path) as fh:
        next(fh)
        for line in fh:
            ddi_id, method, majority_confirmed, mean_confirmed = line.rstrip("\n").split("\t")
            conn.execute(
                "INSERT OR REPLACE INTO ddi_interaction_confirmed "
                "(ddi_id, method, majority_confirmed, mean_confirmed) VALUES (?, ?, ?, ?)",
                (int(ddi_id), method, int(majority_confirmed), int(mean_confirmed)),
            )
            n += 1
    return n


def main():
    args = parse_args()
    shutil.copy(args.db_in, args.db_out)
    conn = sqlite3.connect(args.db_out)

    total_scores = sum(apply_scores(conn, p) for p in args.scores)
    total_confirmed = sum(apply_confirmed(conn, p) for p in args.confirmed)
    conn.commit()

    # Low-contact pairs are skipped by design, so only a method with eligible
    # pairs and not one score is worth a warning (e.g. T_DB=0 in score_ddi.py).
    # The method is in the file name MERGE_SCORE_SHARDS gives: <method>.scores.tsv.
    requested = [os.path.basename(p)[:-len(".scores.tsv")] for p in args.scores
                 if p.endswith(".scores.tsv")]
    counts = {m: (e, s) for m, e, s in conn.execute(
        "SELECT method, COUNT(*), SUM(scored) FROM ("
        "  SELECT method, MAX(z_score IS NOT NULL) AS scored FROM ddi_split_membership"
        "  WHERE has_structure = 1 GROUP BY method, ddi_id, instance_id_a, instance_id_b"
        ") GROUP BY method")}
    for method in requested:
        eligible, scored = counts.get(method, (0, 0))
        print(f"[merge_scores] method={method}: {scored}/{eligible} instance pairs "
              f"with structure scored", flush=True)
        if eligible and not scored:
            print(f"[merge_scores] WARNING: method={method} has {eligible} instance pairs "
                  f"with structure but none scored", file=sys.stderr, flush=True)

    conn.close()
    print(f"[merge_scores] done -> {total_scores} z_scores merged, {total_confirmed} confirmed rows merged", flush=True)

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
