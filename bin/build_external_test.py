#!/usr/bin/env python3
"""Build the external test set: instance pairs for every externally-sourced DDI.

The external sources (``single_domain_ppi``, ``PPIDM*``, ``negatome``) are the
held-out evaluation set.  By the insert rules in ``ddi_db_utils`` a pair that any
protected source claimed was dropped at insert time, so every DDI this script
sees is guaranteed never to have appeared in *any* split of *any* negative set --
strictly stronger than "not in train/validation", and the reason the external
test set needs no ILP, no protein-universe accounting and no leakage analysis of
its own.

For each such DDI ``(A, B)`` it samples up to ``--target`` instance pairs, one
instance of A and one of B, whose **parent proteins are distinct** -- a pair of
domains on the same protein is not evidence of an interaction.  Sampling is
deterministic and independent of iteration order: the RNG for a DDI is seeded
from ``"{seed}:{pfam_a}:{pfam_b}"``, so adding or removing other DDIs never moves
this one's examples.

``--pool-factor`` bounds the work: at most ``target * pool_factor`` instances are
drawn per family before pairs are formed, which caps the cross product at
``(target * pool_factor)^2`` for families with thousands of instances.

The result is written as the ``test`` split of every method named by
``--method``, i.e. both ``external_test`` and ``external_test_hcni`` -- their test
sets are identical by design, the two differ only in their train/validation
negatives.  A DDI that yields no valid pair is dropped and reported.
"""

import argparse
import csv
import random
import sqlite3
import sys
from collections import Counter, defaultdict

from ddi_db_utils import is_protected

DROPPED_COLUMNS = ["pfam_a", "pfam_b", "negative", "source", "reason"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True, help="domainsplit SQLite (modified in place)")
    p.add_argument("--method", required=True, action="append",
                   help="method directory to write this test set under, repeatable")
    p.add_argument("--split", default="test")
    p.add_argument("--target", type=int, required=True, help="instance pairs per DDI (ddi_examples_target)")
    p.add_argument("--pool-factor", type=int, default=1, help="instances drawn per family = target * this")
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--dropped-out", required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def external_ddis(conn):
    """``[(ddi_id, pfam_a, pfam_b, domain_id_a, domain_id_b, negative, source)]`` for non-protected DDIs."""
    rows = []
    for ddi_id, id_a, id_b, pfam_a, pfam_b, negative, source in conn.execute(
        "SELECT ddi.id, ddi.domain_id_a, ddi.domain_id_b, da.pfam_id, db.pfam_id, ddi.negative, ddi.source "
        "FROM domain_domain_interaction AS ddi "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b "
        "ORDER BY ddi.id"
    ):
        if not is_protected(source):
            rows.append((ddi_id, pfam_a, pfam_b, id_a, id_b, negative, source))
    return rows


def instances_by_domain(conn):
    """``{domain_id: [(instance_id, protein_id), ...]}``, sorted for determinism."""
    by_domain = defaultdict(list)
    for domain_id, instance_id, protein_id in conn.execute(
        "SELECT domain_id, instance_id, protein_id FROM domain_protein_map "
        "WHERE instance_id IS NOT NULL ORDER BY instance_id"
    ):
        by_domain[domain_id].append((instance_id, protein_id))
    return by_domain


def draw_pool(rng, instances, size):
    """``size`` instances drawn without replacement, or all of them if there are fewer."""
    if len(instances) <= size:
        return list(instances)
    return rng.sample(instances, size)


def sample_pairs(rng, inst_a, inst_b, target, pool_size, same_family):
    """Up to ``target`` unordered instance pairs with distinct parent proteins."""
    pool_a = draw_pool(rng, inst_a, pool_size)
    pool_b = pool_a if same_family else draw_pool(rng, inst_b, pool_size)

    candidates = set()
    for id_a, parent_a in pool_a:
        for id_b, parent_b in pool_b:
            if parent_a == parent_b:
                continue
            first, second = sorted((id_a, id_b))
            candidates.add((first, second))
    if not candidates:
        return []
    ordered = sorted(candidates)
    rng.shuffle(ordered)
    return ordered[:target]


def main():
    args = parse_args()

    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    ddis = external_ddis(conn)
    by_domain = instances_by_domain(conn)
    pool_size = max(1, args.target * max(1, args.pool_factor))

    membership, dropped, stats = [], [], Counter()
    for ddi_id, pfam_a, pfam_b, id_a, id_b, negative, source in ddis:
        stats["ddis"] += 1
        inst_a, inst_b = by_domain.get(id_a, []), by_domain.get(id_b, [])
        if not inst_a or not inst_b:
            dropped.append((pfam_a, pfam_b, negative, source, "no_instances"))
            stats["dropped_no_instances"] += 1
            continue
        rng = random.Random(f"{args.seed}:{pfam_a}:{pfam_b}")
        pairs = sample_pairs(rng, inst_a, inst_b, args.target, pool_size, same_family=(id_a == id_b))
        if not pairs:
            dropped.append((pfam_a, pfam_b, negative, source, "no_distinct_parent_pair"))
            stats["dropped_no_distinct_parent"] += 1
            continue
        if len(pairs) < args.target:
            stats["partial"] += 1
        stats["pairs"] += len(pairs)
        for method in args.method:
            for first, second in pairs:
                membership.append((ddi_id, method, args.split, first, second))

    conn.executemany(
        "INSERT OR IGNORE INTO ddi_split_membership"
        "(ddi_id, method, split, instance_id_a, instance_id_b) VALUES (?, ?, ?, ?, ?)",
        membership,
    )
    conn.commit()

    with open(args.dropped_out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(DROPPED_COLUMNS)
        writer.writerows(sorted(dropped))

    for key in ("ddis", "pairs", "partial", "dropped_no_instances", "dropped_no_distinct_parent"):
        print(f"[external_test] {key} = {stats[key]}", flush=True)
    print(f"[external_test] methods = {','.join(args.method)} split = {args.split}", flush=True)
    conn.close()

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
