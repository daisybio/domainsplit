#!/usr/bin/env python3
"""Delete DDIs that no domain instance can represent, and the rows they stranded.

A Pfam family can end up with zero instances: ``--tier human_only`` drops
families whose strata are all non-human, and a dead or mistyped accession never
resolves at all.  A DDI touching such a family cannot be turned into an instance
pair, so it cannot be split, cannot be embedded, and cannot appear in the
external test set -- it would sit in ``domainsplit.sqlite3`` belonging to
nothing.

One module owns the invariant "every DDI in the master DB is representable".
That is what lets the inserters stay simple: they keep bulk-creating missing
``domain`` rows without having to know which families have instances, and this
step cleans up afterwards.

Deleted, in order:

1. every DDI where either family has no ``domain_protein_map`` row, whatever its
   source -- 3did included, since the orphan 3did DDIs that ``human_only``
   strands are exactly this case;
2. every ``domain`` row no surviving DDI references (this cascades to its
   ``domain_protein_map`` rows);
3. every ``protein`` row no surviving mapping references, so ENRICH does not
   fetch sequences and embeddings for proteins nothing will ever ask about.

``ddi_split_membership`` rows follow their DDI through ``ON DELETE CASCADE``.
"""

import argparse
import csv
import sqlite3
import sys
from collections import Counter

from ddi_db_utils import split_sources

PRUNE_COLUMNS = ["pfam_a", "pfam_b", "negative", "source", "reason"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True, help="domainsplit SQLite (modified in place)")
    p.add_argument("--report-out", required=True, help="pruned_ddis.tsv")
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def unrepresented_ddis(conn):
    """``[(ddi_id, pfam_a, pfam_b, negative, source, reason)]`` for DDIs with a family that has no instance."""
    with_instances = {
        row[0] for row in conn.execute("SELECT DISTINCT domain_id FROM domain_protein_map")
    }
    doomed = []
    for ddi_id, id_a, id_b, pfam_a, pfam_b, negative, source in conn.execute(
        "SELECT ddi.id, ddi.domain_id_a, ddi.domain_id_b, da.pfam_id, db.pfam_id, ddi.negative, ddi.source "
        "FROM domain_domain_interaction AS ddi "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b"
    ):
        missing = [p for p, i in ((pfam_a, id_a), (pfam_b, id_b)) if i not in with_instances]
        if missing:
            doomed.append((ddi_id, pfam_a, pfam_b, negative, source, "no_instances:" + ",".join(missing)))
    return doomed


def main():
    args = parse_args()

    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    n_ddis_before = conn.execute("SELECT COUNT(*) FROM domain_domain_interaction").fetchone()[0]
    doomed = unrepresented_ddis(conn)

    per_source = Counter()
    for _, _, _, _, source, _ in doomed:
        for s in split_sources(source):
            per_source[s] += 1

    conn.executemany(
        "DELETE FROM domain_domain_interaction WHERE id = ?", [(row[0],) for row in doomed]
    )
    # Domains no surviving DDI references; the delete cascades to their instances.
    n_domains = conn.execute(
        "DELETE FROM domain WHERE id NOT IN ("
        "  SELECT domain_id_a FROM domain_domain_interaction"
        "  UNION SELECT domain_id_b FROM domain_domain_interaction)"
    ).rowcount
    n_proteins = conn.execute(
        "DELETE FROM protein WHERE id NOT IN (SELECT protein_id FROM domain_protein_map)"
    ).rowcount
    conn.commit()

    with open(args.report_out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(PRUNE_COLUMNS)
        writer.writerows(sorted(row[1:] for row in doomed))

    n_ddis_after = conn.execute("SELECT COUNT(*) FROM domain_domain_interaction").fetchone()[0]
    print(f"[prune] ddis_before = {n_ddis_before}", flush=True)
    print(f"[prune] ddis_pruned = {len(doomed)}", flush=True)
    print(f"[prune] ddis_after = {n_ddis_after}", flush=True)
    print(f"[prune] domains_deleted = {n_domains}", flush=True)
    print(f"[prune] proteins_deleted = {n_proteins}", flush=True)
    for source in sorted(per_source):
        print(f"[prune] pruned_source_{source} = {per_source[source]}", flush=True)
    conn.close()

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
