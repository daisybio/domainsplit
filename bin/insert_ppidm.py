#!/usr/bin/env python3
"""Insert PPIDM predicted positive DDIs, keeping the class as the source.

Input ``predicted_ddi_ppi.tsv`` columns: ``domain_1  domain_2  class`` where each
domain token looks like ``10114/PF00069`` (the Pfam accession follows the slash).
Rows are inserted as ``negative=0, source='PPIDM_<Class>'`` for the requested
classes.  Classes are processed Gold -> Silver -> Bronze so that, on a duplicate
domain pair, the highest-confidence class wins (``INSERT OR IGNORE``).
"""

import argparse
import sqlite3
import sys

from ddi_db_utils import count_source, ensure_domains, insert_ddis

# highest confidence first
CLASS_ORDER = ["Gold", "Silver", "Bronze"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True)
    p.add_argument("--ppidm", required=True)
    p.add_argument("--classes", required=True,
                   help="comma-separated classes to include, e.g. 'Bronze,Silver,Gold'")
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def extract_pfam(token):
    """``10114/PF00069`` or ``PF00069.3`` -> ``PF00069`` (or None)."""
    pf = token.split("/")[-1].split(".")[0].strip()
    return pf if pf.startswith("PF") else None


def main():
    args = parse_args()

    allowed = {c.strip().capitalize() for c in args.classes.split(",") if c.strip()}
    classes = [c for c in CLASS_ORDER if c in allowed]
    print(f"[ppidm] including classes: {classes}", flush=True)

    # class -> list of (pfam_a, pfam_b)
    pairs_by_class = {c: [] for c in classes}
    n_rows = n_bad = 0
    with open(args.ppidm) as fh:
        for i, line in enumerate(fh):
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            if i == 0 and cols[2].strip().lower() == "class":
                continue  # header
            cls = cols[2].strip().capitalize()
            if cls not in pairs_by_class:
                continue
            pfam_a = extract_pfam(cols[0])
            pfam_b = extract_pfam(cols[1])
            if pfam_a is None or pfam_b is None:
                n_bad += 1
                continue
            pairs_by_class[cls].append((pfam_a, pfam_b))
            n_rows += 1

    print(f"[ppidm] parsed {n_rows} pairs ({n_bad} unparseable)", flush=True)

    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    all_pfams = {p for pairs in pairs_by_class.values() for pair in pairs for p in pair}
    ensure_domains(conn, all_pfams)

    for cls in classes:  # Gold first
        source = f"PPIDM_{cls}"
        insert_ddis(conn, pairs_by_class[cls], negative=False, source=source)
        conn.commit()
        print(f"[ppidm] n_ddis_source_{source} = {count_source(conn, source)}", flush=True)

    conn.close()

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
