#!/usr/bin/env python3
"""
merge_source_scores.py
-----------------------
Phase 3 of DDI scoring: recombine the per-source scored DBs for one
domain-splitting method back into a single DB.

score_ddi.py now runs once per (split x source) and only ever writes its
own source's 2 columns (interaction_confirmed_majority_<source>,
interaction_confirmed_mean_<source>). This script takes the AF3-scored DB
and the RF-scored DB for the SAME split and copies RF's 2 columns onto a
copy of the AF3 DB, matched by domain_domain_interaction.id, so downstream
consumers see all 4 columns on one DB per split.
"""

import argparse
import shutil
import sqlite3
import sys


SOURCES = ("AF3", "RF")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--db_af3", required=True, help="dbscored.sqlite3 from the AF3 SCORE_DDI run")
    p.add_argument("--db_rf", required=True, help="dbscored.sqlite3 from the RF SCORE_DDI run")
    p.add_argument("--db_out", required=True, help="Merged output DB path")
    p.add_argument("--versions", required=True)
    p.add_argument("--process_name", required=True)
    return p.parse_args()


def get_source_columns(conn, source):
    cols = [f"interaction_confirmed_majority_{source}", f"interaction_confirmed_mean_{source}"]
    existing = {row[1] for row in conn.execute("PRAGMA table_info(domain_domain_interaction)")}
    present = all(c in existing for c in cols)
    return cols, present


def merge(db_af3, db_rf, db_out):
    # Base DB = AF3 output (already has the AF3 columns); copy RF's columns in.
    shutil.copy(db_af3, db_out)
    conn = sqlite3.connect(db_out)
    conn.execute("PRAGMA foreign_keys = ON")

    af3_cols, af3_present = get_source_columns(conn, "AF3")
    if af3_present:
        print(f"[merge_sources] AF3 columns already present in {db_out}, skipping AF3 merge", flush=True)
    else:
        print(f"[merge_sources] WARNING: AF3 columns not present in {db_out}, no AF3 predictions available", flush=True)

    rf_cols = [f"interaction_confirmed_majority_RF", f"interaction_confirmed_mean_RF"]
    for col in rf_cols:
        conn.execute(f"ALTER TABLE domain_domain_interaction ADD COLUMN {col} INTEGER DEFAULT 0;")

    conn.execute("ATTACH DATABASE ? AS rf_db", (db_rf,))
    rf_cols = ["interaction_confirmed_majority_RF", "interaction_confirmed_mean_RF"]
    rf_check_cols = {row[1] for row in conn.execute("PRAGMA rf_db.table_info(domain_domain_interaction)")}
    rf_present = all(c in rf_check_cols for c in rf_cols)

    n_updated = 0
    if rf_present:
        for col in rf_cols:
            conn.execute(f"ALTER TABLE domain_domain_interaction ADD COLUMN {col} INTEGER;")
        conn.execute(f"""
            UPDATE domain_domain_interaction
            SET interaction_confirmed_majority_RF = (
                    SELECT rf.interaction_confirmed_majority_RF
                    FROM rf_db.domain_domain_interaction rf
                    WHERE rf.id = domain_domain_interaction.id
                ),
                interaction_confirmed_mean_RF = (
                    SELECT rf.interaction_confirmed_mean_RF
                    FROM rf_db.domain_domain_interaction rf
                    WHERE rf.id = domain_domain_interaction.id
                )
            WHERE id IN (SELECT id FROM rf_db.domain_domain_interaction)
        """)
        n_updated = conn.total_changes
    else:
        print(f"[merge_sources] RF columns absent from {db_rf} "
              f"(no RF predictions for this split) — leaving them out entirely", flush=True)

    conn.commit()
    conn.execute("DETACH DATABASE rf_db")
    conn.close()
    
    print(f"[merge_sources] merged RF columns into {db_out} "
          f"({n_updated} rows touched)", flush=True)


def write_versions(path, process_name):
    with open(path, "w") as f:
        f.write(f'"{process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")


def main():
    args = parse_args()
    merge(args.db_af3, args.db_rf, args.db_out)
    write_versions(args.versions, args.process_name)
    print("[merge_sources] done", flush=True)


if __name__ == "__main__":
    main()