#!/usr/bin/env python3
"""
merge_source_scores.py
-----------------------
Phase 3 of DDI scoring: recombine the per-source scored DBs for one split
back into a single DB. Handles any number of configured sources (1..N).

score_ddi.py runs once per (split x source) and only ever writes its own
source's 2 columns (interaction_confirmed_majority_<source>,
interaction_confirmed_mean_<source>) -- and only if that source had any
predicted structures at all for this split. This script takes the scored
DB for each configured source and copies each source's columns onto a
single base DB, matched by domain_domain_interaction.id. Sources with no
columns present in their DB (no predictions) are skipped entirely -- no
empty columns are created for them.
"""

import argparse
import shutil
import sqlite3
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--db", action="append", required=True, metavar="SOURCE=PATH",
        help="One scored DB per configured source, e.g. --db AF3=af3.sqlite3 "
             "--db RF=rf.sqlite3. Repeatable; at least one required."
    )
    p.add_argument("--db_out", required=True, help="Merged output DB path")
    p.add_argument("--versions", required=True)
    p.add_argument("--process_name", required=True)
    return p.parse_args()


def parse_source_dbs(raw):
    parsed = []
    for entry in raw:
        if "=" not in entry:
            raise ValueError(f"--db entries must be SOURCE=PATH, got: {entry}")
        source, path = entry.split("=", 1)
        parsed.append((source, path))
    return parsed


def get_source_columns(conn, source, schema_prefix=None):
    cols = [f"interaction_confirmed_majority_{source}", f"interaction_confirmed_mean_{source}"]
    table = f"{schema_prefix}.domain_domain_interaction" if schema_prefix else "domain_domain_interaction"
    existing = {row[1] for row in conn.execute(f"PRAGMA table_info({table})")}
    return cols, all(c in existing for c in cols)


def merge(source_dbs, db_out):
    if not source_dbs:
        raise ValueError("At least one --db SOURCE=PATH is required")

    base_source, base_path = source_dbs[0]
    shutil.copy(base_path, db_out)

    if len(source_dbs) == 1:
        print(f"[merge_sources] only one source ({base_source}) -- no merge needed", flush=True)
        return

    
    conn = sqlite3.connect(db_out)
    conn.execute("PRAGMA foreign_keys = ON")

    base_cols, base_present = get_source_columns(conn, base_source)
    if base_present:
        print(f"[merge_sources] base DB ({base_source}) already has columns: {base_cols}", flush=True)
    else:
        print(f"[merge_sources] WARNING: base DB ({base_source}) has no columns -- "
              f"no {base_source} predictions for this split", flush=True)

    for i, (source, path) in enumerate(source_dbs[1:], start=1):
        schema = f"src_{i}"
        conn.execute("ATTACH DATABASE ? AS ?", (path, schema))
        cols, present = get_source_columns(conn, source, schema_prefix=schema)

        if not present:
            print(f"[merge_sources] {source} columns absent from {path} "
                  f"(no {source} predictions for this split) -- leaving them out entirely", flush=True)
            conn.execute(f"DETACH DATABASE {schema}")
            continue

        for col in cols:
            conn.execute(f"ALTER TABLE domain_domain_interaction ADD COLUMN {col} INTEGER;")

        conn.execute(f"""
            UPDATE domain_domain_interaction
            SET {cols[0]} = (
                    SELECT src.{cols[0]}
                    FROM {schema}.domain_domain_interaction src
                    WHERE src.id = domain_domain_interaction.id
                ),
                {cols[1]} = (
                    SELECT src.{cols[1]}
                    FROM {schema}.domain_domain_interaction src
                    WHERE src.id = domain_domain_interaction.id
                )
            WHERE id IN (SELECT id FROM {schema}.domain_domain_interaction)
        """)
        n_updated = conn.total_changes
        conn.execute(f"DETACH DATABASE {schema}")
        print(f"[merge_sources] merged {source} columns into {db_out} ({n_updated} rows touched)", flush=True)

    conn.commit()
    conn.close()


def write_versions(path, process_name):
    with open(path, "w") as f:
        f.write(f'"{process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")


def main():
    args = parse_args()
    source_dbs = parse_source_dbs(args.db)
    merge(source_dbs, args.db_out)
    write_versions(args.versions, args.process_name)
    print("[merge_sources] done", flush=True)


if __name__ == "__main__":
    main()