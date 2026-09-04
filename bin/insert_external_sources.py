#!/usr/bin/env python3
"""Insert every external DDI source into the domainsplit SQLite, in one pass.

Consumes the normalized TSVs written by the ``PARSE_*`` modules and applies the
merge/drop rules in ``ddi_db_utils.merge_external_ddis``:

  * a pair already held by a protected source (``3did``, ``sampled_negative``)
    is dropped silently -- that is what makes the external test set strictly
    unseen;
  * external sources agreeing on ``negative`` merge into one row whose ``source``
    is their comma-joined list;
  * external sources disagreeing on ``negative`` drop the pair entirely, and any
    existing external row for it is deleted.

Every dropped contribution is written to the conflict report, one row per
(pair, source), so the report can be grepped by source.  The whole batch is
resolved before anything is written, so the outcome does not depend on the order
the TSVs are supplied in.
"""

import argparse
import csv
import sqlite3
import sys

from ddi_db_utils import CONFLICT_COLUMNS, count_ddis, count_source, merge_external_ddis
from external_ddi_tsv import read_external_tsv


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True, help="domainsplit SQLite (modified in place)")
    p.add_argument("--ddis", required=True, nargs="+",
                   help="normalized external-DDI TSVs (pfam_a, pfam_b, negative, source)")
    p.add_argument("--conflicts-out", required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    rows = []
    for path in args.ddis:
        n_before = len(rows)
        rows.extend(read_external_tsv(path))
        print(f"[external] {path}: {len(rows) - n_before} rows", flush=True)
    print(f"[external] {len(rows)} rows over {len(args.ddis)} sources", flush=True)

    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    stats, conflicts = merge_external_ddis(conn, rows)
    conn.commit()

    with open(args.conflicts_out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(CONFLICT_COLUMNS)
        writer.writerows(conflicts)

    for key in ("pairs_offered", "inserted", "merged", "dropped_protected",
                "conflict_pairs", "sources_dropped"):
        print(f"[external] {key} = {stats[key]}", flush=True)
    for source in sorted({row[3] for row in rows}):
        print(f"[external] n_ddis_source_{source} = {count_source(conn, source)}",
              flush=True)
    print(f"[external] n_ddis_total = {count_ddis(conn)}", flush=True)
    conn.close()

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
