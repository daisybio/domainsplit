#!/usr/bin/env python3
"""Insert ppi-splitting's sampled negative DDIs into the domainsplit SQLite.

Every family pair carrying ``label=0`` in any family-level split CSV, of any
split, of any negative set, becomes one row with ``negative=1`` and
``source='sampled_negative'``.  The pair is stored once no matter how many splits
or negative sets produced it; which split it belongs to is recorded separately by
INGEST_SPLIT_MEMBERSHIP.

``sampled_negative`` is a *protected* source (see ``ddi_db_utils``): it cannot
collide with ``3did`` by construction -- ppi-splitting samples negatives from
pairs that are not positives -- and it owns its pairs against every external
source.  That only holds if this step runs **before** INSERT_EXTERNAL_SOURCES, so
the script counts pairs it had to skip because an external row already held them
and fails if there are any: a non-zero count means the insertion order
``3did -> sampled_negative -> external`` was broken, and the external test set
would no longer be strictly unseen.
"""

import argparse
import sqlite3
import sys
from collections import Counter

from ddi_db_utils import canonical_pair, count_source, ensure_domains, insert_ddis, is_protected
from split_io import parse_split_spec, read_family_csv

SOURCE = "sampled_negative"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True, help="domainsplit SQLite (modified in place)")
    p.add_argument("--split", required=True, action="append", metavar="METHOD:SPLIT:PATH",
                   help="family-level split CSV, repeatable")
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def collect_negatives(specs):
    """``(pairs, stats)`` -- the canonicalised set of ``label=0`` family pairs."""
    pairs, stats = set(), Counter()
    for spec in specs:
        method, split, path = parse_split_spec(spec)
        n = 0
        for pfam_a, pfam_b, label in read_family_csv(path):
            stats["rows"] += 1
            if label == 0:
                pairs.add(canonical_pair(pfam_a, pfam_b))
                n += 1
        stats["negative_rows"] += n
        print(f"[negatives] {method}/{split}: {n} negative rows", flush=True)
    return pairs, stats


def existing_owners(conn, pairs):
    """``{pair: source}`` for those ``pairs`` that already have a DDI row."""
    out = {}
    for pfam_a, pfam_b, source in conn.execute(
        "SELECT da.pfam_id, db.pfam_id, ddi.source "
        "FROM domain_domain_interaction AS ddi "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b"
    ):
        pair = canonical_pair(pfam_a, pfam_b)
        if pair in pairs:
            out[pair] = source
    return out


def main():
    args = parse_args()

    pairs, stats = collect_negatives(args.split)
    print(f"[negatives] {len(pairs)} distinct pairs from {stats['negative_rows']} rows", flush=True)

    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    prior = existing_owners(conn, pairs)
    protected_prior = {p: s for p, s in prior.items() if is_protected(s)}
    external_prior = {p: s for p, s in prior.items() if not is_protected(s)}

    ensure_domains(conn, {pfam for pair in pairs for pfam in pair})
    offered = insert_ddis(conn, sorted(pairs), negative=True, source=SOURCE)
    conn.commit()

    print(f"[negatives] offered = {offered}", flush=True)
    print(f"[negatives] already_held_by_protected = {len(protected_prior)}", flush=True)
    print(f"[negatives] n_ddis_sampled_negative = {count_source(conn, SOURCE)}", flush=True)
    conn.close()

    if external_prior:
        sample = ", ".join(f"{a}/{b} ({s})" for (a, b), s in sorted(external_prior.items())[:10])
        raise SystemExit(
            f"ERROR: {len(external_prior)} sampled-negative pairs are already held by an external "
            f"source: {sample}" + ("..." if len(external_prior) > 10 else "") + "\n"
            "INSERT_EXTERNAL_SOURCES must run after this step -- the insertion order "
            "3did -> sampled_negative -> external is what makes the external test set unseen."
        )

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
