#!/usr/bin/env python3
"""Merge every method's scores.tsv/confirmed.tsv (score_ddi.py) into one copy
of the struct-enriched master. Each method's files touch only its own
ddi_split_membership rows and ddi_interaction_confirmed rows, so the merges
are disjoint and order-independent."""
import argparse as ap
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
            ddi_id, method, split, ia, ib, z_score = line.rstrip("\n").split("\t")
            conn.execute(
                "UPDATE ddi_split_membership SET z_score = ? "
                "WHERE ddi_id = ? AND method = ? AND split = ? "
                "AND instance_id_a = ? AND instance_id_b = ?",
                (float(z_score), int(ddi_id), method, split, ia, ib),
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
                (ddi_id, method, int(majority_confirmed), int(mean_confirmed)),
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

    unscored = conn.execute("SELECT COUNT(*) FROM ddi_split_membership WHERE z_score IS NULL").fetchone()[0]
    scored = conn.execute("SELECT COUNT(*) FROM ddi_split_membership WHERE z_score IS NOT NULL").fetchone()[0]
    if unscored:
        print(f"[merge_scores] WARNING: {unscored} ddi_split_membership rows still unscored", flush=True)

    print(f"[merge_scores] {scored} ddi_split_membership rows scored", flush=True)
    
    conn.close()
    print(f"[merge_scores] done -> {total_scores} z_scores merged, {total_confirmed} confirmed rows merged", flush=True)

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()