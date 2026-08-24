#!/usr/bin/env python3
"""Reduce a copy of the enriched master DB to one ``(method, split)``.

``ddi_split_membership`` already records which DDIs belong to which split, so
this is a pure SQL filter: keep that split's DDIs, then delete everything no
surviving DDI reaches.  One module therefore replaces all three of the old
splitter modules -- the leakage-reduction logic lives in ppi-splitting, and what
is left here is bookkeeping.

Deletion order matters, and each step relies on the previous one having run:

1. DDIs not in this split;
2. ``domain`` rows no surviving DDI references -- cascades to their
   ``domain_protein_map`` and ``domain_go_terms`` rows;
3. ``protein`` rows no surviving mapping references -- cascades to their
   ``protein_go_terms`` and ``protein_protein_interaction`` rows;
4. ``ddi_split_membership`` rows for every other ``(method, split)``, so each
   published DB describes only itself.

Instances of a surviving family are kept whole, not narrowed to the pairs this
split happens to use: the split is a partition of *families*, so an instance of a
kept family belongs to this split by construction, and keeping them lets a
consumer form its own pairs.
"""

import argparse
import sqlite3
import sys


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True, help="copy of the enriched master DB, modified in place")
    p.add_argument("--method", required=True)
    p.add_argument("--split", required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def counts(conn):
    tables = ["domain_domain_interaction", "domain", "protein", "domain_protein_map",
              "domain_go_terms", "protein_go_terms", "protein_protein_interaction"]
    return {t: conn.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0] for t in tables}


def main():
    args = parse_args()

    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    before = counts(conn)
    kept = conn.execute(
        "SELECT COUNT(DISTINCT ddi_id) FROM ddi_split_membership WHERE method = ? AND split = ?",
        (args.method, args.split),
    ).fetchone()[0]
    if kept == 0:
        print(
            f"[subset] WARNING: no ddi_split_membership rows for {args.method}/{args.split}; "
            "the emitted database will be empty.",
            flush=True,
        )

    conn.execute(
        "DELETE FROM domain_domain_interaction WHERE id NOT IN ("
        "  SELECT ddi_id FROM ddi_split_membership WHERE method = ? AND split = ?)",
        (args.method, args.split),
    )
    conn.execute(
        "DELETE FROM domain WHERE id NOT IN ("
        "  SELECT domain_id_a FROM domain_domain_interaction"
        "  UNION SELECT domain_id_b FROM domain_domain_interaction)"
    )
    conn.execute("DELETE FROM protein WHERE id NOT IN (SELECT protein_id FROM domain_protein_map)")
    conn.execute(
        "DELETE FROM ddi_split_membership WHERE method != ? OR split != ?", (args.method, args.split)
    )
    conn.commit()

    after = counts(conn)
    print(f"[subset] {args.method}/{args.split}: {kept} DDIs kept", flush=True)
    for table in before:
        print(f"[subset] {table} = {after[table]} (was {before[table]})", flush=True)

    conn.isolation_level = None
    conn.execute("VACUUM")
    conn.close()

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
