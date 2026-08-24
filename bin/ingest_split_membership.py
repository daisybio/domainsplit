#!/usr/bin/env python3
"""Record which DDIs and which instance pairs belong to which ``(method, split)``.

Reads ppi-splitting's instance-level split CSVs (``protein1,protein2,label,
family1,family2``, where ``protein1``/``protein2`` are instance ids) and fills
``ddi_split_membership``.  That table exists for one reason: it turns
SUBSET_SPLIT_DB into a pure SQL filter, so a single module can replace all three
of the old splitter modules.

Both labels are recorded.  A ``label=1`` pair is a 3did DDI, a ``label=0`` pair
is one INGEST_SAMPLED_NEGATIVES has just inserted; either way the family pair
already has a row in ``domain_domain_interaction`` by the time this runs, so a
row that fails to resolve is a wiring error, not data.  Unresolved rows are
counted and written to a report rather than dropped silently.

The ``val`` -> ``validation`` rename happens here (see ``split_io``).
"""

import argparse
import csv
import sqlite3
import sys
from collections import Counter

from ddi_db_utils import canonical_pair
from split_io import parse_split_spec, read_instance_csv

UNRESOLVED_COLUMNS = ["method", "split", "family_a", "family_b", "instance_a", "instance_b", "label", "reason"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True, help="domainsplit SQLite (modified in place)")
    p.add_argument("--split", required=True, action="append", metavar="METHOD:SPLIT:PATH",
                   help="instance-level split CSV (*_instances.csv), repeatable")
    p.add_argument("--unresolved-out", required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def ddi_ids_by_pair(conn):
    """``{(pfam_a, pfam_b): ddi_id}`` keyed by canonicalised pair."""
    return {
        canonical_pair(a, b): ddi_id
        for ddi_id, a, b in conn.execute(
            "SELECT ddi.id, da.pfam_id, db.pfam_id "
            "FROM domain_domain_interaction AS ddi "
            "JOIN domain AS da ON da.id = ddi.domain_id_a "
            "JOIN domain AS db ON db.id = ddi.domain_id_b"
        )
    }


def known_instances(conn):
    return {row[0] for row in conn.execute(
        "SELECT instance_id FROM domain_protein_map WHERE instance_id IS NOT NULL"
    )}


def collect(specs, ddi_id_of, instances):
    """``(membership_rows, unresolved_rows, stats)``.

    Instance pairs are canonicalised alongside their families so a pair written
    as ``(B, A)`` in one split and ``(A, B)`` in another produces one row.
    """
    rows, unresolved, stats = set(), [], Counter()
    for spec in specs:
        method, split, path = parse_split_spec(spec)
        n = 0
        for inst_a, inst_b, label, fam_a, fam_b in read_instance_csv(path):
            stats["rows"] += 1
            ddi_id = ddi_id_of.get(canonical_pair(fam_a, fam_b))
            if ddi_id is None:
                unresolved.append((method, split, fam_a, fam_b, inst_a, inst_b, label, "no_ddi_row"))
                stats["no_ddi_row"] += 1
                continue
            missing = [i for i in (inst_a, inst_b) if i not in instances]
            if missing:
                unresolved.append((method, split, fam_a, fam_b, inst_a, inst_b, label, "unknown_instance"))
                stats["unknown_instance"] += 1
                continue
            first, second = sorted((inst_a, inst_b))
            rows.add((ddi_id, method, split, first, second))
            n += 1
        stats["recorded"] += n
        print(f"[membership] {method}/{split}: {n} instance pairs", flush=True)
    return sorted(rows), unresolved, stats


def main():
    args = parse_args()

    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    rows, unresolved, stats = collect(args.split, ddi_ids_by_pair(conn), known_instances(conn))

    conn.executemany(
        "INSERT OR IGNORE INTO ddi_split_membership"
        "(ddi_id, method, split, instance_id_a, instance_id_b) VALUES (?, ?, ?, ?, ?)",
        rows,
    )
    conn.commit()

    with open(args.unresolved_out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(UNRESOLVED_COLUMNS)
        writer.writerows(unresolved)

    for key in ("rows", "recorded", "no_ddi_row", "unknown_instance"):
        print(f"[membership] {key} = {stats[key]}", flush=True)
    for method, split, n, n_ddis in conn.execute(
        "SELECT method, split, COUNT(*), COUNT(DISTINCT ddi_id) "
        "FROM ddi_split_membership GROUP BY method, split ORDER BY method, split"
    ):
        print(f"[membership] {method}/{split} = {n} pairs over {n_ddis} DDIs", flush=True)
    conn.close()

    if stats["no_ddi_row"] or stats["unknown_instance"]:
        print(
            f"[membership] WARNING: {stats['no_ddi_row'] + stats['unknown_instance']} rows did not "
            f"resolve; see {args.unresolved_out}. Every split row should already have a DDI row "
            "(3did for label=1, sampled_negative for label=0) and both instances in "
            "domain_protein_map, so this points at a wiring problem, not at the data.",
            flush=True,
        )

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
