#!/usr/bin/env python3
"""Insert Negatome negative DDIs into the domainsplit SQLite.

Negatome ``combined_pfam.txt`` lists whitespace-separated Pfam pairs that do not
interact.  Treated like every other source: bulk-create missing domains, then
insert as ``negative=1, source='negatome'``.
"""

import argparse
import sqlite3
import sys

from ddi_db_utils import count_source, ensure_domains, insert_ddis


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True)
    p.add_argument("--negatome", required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def iter_negatome_pairs(path):
    with open(path) as f:
        for line in f:
            tokens = line.split()
            if len(tokens) < 2:
                continue
            yield tokens[0], tokens[1]


def main():
    args = parse_args()

    pairs = list(iter_negatome_pairs(args.negatome))
    print(f"[negatome] read {len(pairs)} pairs", flush=True)

    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    pfams = {p for pair in pairs for p in pair}
    ensure_domains(conn, pfams)
    insert_ddis(conn, pairs, negative=True, source="negatome")
    conn.commit()
    print(f"[negatome] n_ddis_source_negatome = {count_source(conn, 'negatome')}", flush=True)
    conn.close()

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
