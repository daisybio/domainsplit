#!/usr/bin/env python3
"""Insert 3did positive DDIs into the domainsplit SQLite.

3did is treated like every other source: read its DDI pairs, bulk-create any
missing ``domain`` rows for the referenced Pfam IDs, then insert the
interactions as ``negative=0, source='3did'``.
"""

import argparse
import sqlite3
import sys

from ddi_db_utils import count_source, ensure_domains, insert_ddis


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True, help="domainsplit SQLite (modified in place)")
    p.add_argument("--sqlite-3did", required=True, help="3did SQLite from DOWNLOAD_3DID_SQLITE")
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def iter_3did_pairs(conn_3did):
    """Yield (pfam_a, pfam_b) Pfam accessions (version stripped) for each 3did DDI."""
    cursor = conn_3did.execute(
        "SELECT d1.Pfam_id, d2.Pfam_id "
        "FROM DDI1, Domain AS d1, Domain AS d2 "
        "WHERE DDI1.domain1 = d1.Name AND DDI1.domain2 = d2.Name"
    )
    for id_1, id_2 in cursor:
        yield id_1.split(".")[0], id_2.split(".")[0]


def main():
    args = parse_args()

    conn_3did = sqlite3.connect(args.sqlite_3did)
    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    pairs = list(iter_3did_pairs(conn_3did))
    conn_3did.close()
    print(f"[3did] read {len(pairs)} DDI pairs", flush=True)

    pfams = {p for pair in pairs for p in pair}
    n_domains = ensure_domains(conn, pfams)
    print(f"[3did] ensured {n_domains} domains", flush=True)

    insert_ddis(conn, pairs, negative=False, source="3did")
    conn.commit()
    print(f"[3did] n_ddis_source_3did = {count_source(conn, '3did')}", flush=True)
    conn.close()

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
