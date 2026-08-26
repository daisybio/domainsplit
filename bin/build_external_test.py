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
import time
from collections import Counter, defaultdict

from ddi_db_utils import is_protected

DROPPED_COLUMNS = ["pfam_a", "pfam_b", "negative", "source", "reason"]

_T0 = time.monotonic()


def log(msg):
    """Progress line with elapsed seconds.

    This process printed nothing at all until its last line, so its first real
    run -- which sat RUNNING for three hours on the cluster and had to be
    cancelled -- left no way to tell the DB copy from the connect from the scan.
    Every stage now announces itself, unbuffered.
    """
    print(f"[external_test] +{time.monotonic() - _T0:7.1f}s {msg}", flush=True)


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
    """Yield ``(ddi_id, pfam_a, pfam_b, domain_id_a, domain_id_b, negative, source)``
    for non-protected DDIs, in ``ddi.id`` order.

    A generator, not a list: at production scale the DDI table is dominated by
    PPIDM (millions of rows), and materialising it only to walk it once cost a
    gigabyte for nothing.  ``is_protected`` has to stay in Python -- ``source``
    can hold a comma-joined list (``"PPIDM,PPIDM_Gold"``) that SQL cannot split.
    """
    for ddi_id, id_a, id_b, pfam_a, pfam_b, negative, source in conn.execute(
        "SELECT ddi.id, ddi.domain_id_a, ddi.domain_id_b, da.pfam_id, db.pfam_id, ddi.negative, ddi.source "
        "FROM domain_domain_interaction AS ddi "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b "
        "ORDER BY ddi.id"
    ):
        if not is_protected(source):
            yield (ddi_id, pfam_a, pfam_b, id_a, id_b, negative, source)


def instances_by_domain(conn):
    """``{domain_id: [(instance_id, protein_id), ...]}``, sorted for determinism.

    Bounded by design rather than by the protein universe: ``FETCH_DOMAIN_META``
    samples at most ``ddi_examples_target * ddi_examples_pool_factor`` instances
    per family, so this is families x that factor, not one row per domain
    occurrence in the proteome.
    """
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

    log(f"connecting to {args.db}")
    conn = sqlite3.connect(args.db)
    # Fail on a lock rather than block on one: an indefinite wait here is
    # indistinguishable from slow work, which is exactly what made the
    # three-hour stall unreadable.
    conn.execute("PRAGMA busy_timeout=120000")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    log("connected, pragmas set")

    by_domain = instances_by_domain(conn)
    log(f"loaded instances for {len(by_domain)} domains "
        f"({sum(len(v) for v in by_domain.values())} instances)")
    pool_size = max(1, args.target * max(1, args.pool_factor))

    # One row per (DDI, method, instance pair): target * len(methods) times the
    # external DDI count, which is millions once PPIDM is in. Flushed in batches
    # so peak memory is the batch, not the whole insert -- and on a separate
    # cursor, because the read cursor above is still streaming the DDI table.
    writer_cur = conn.cursor()
    INSERT = ("INSERT OR IGNORE INTO ddi_split_membership"
              "(ddi_id, method, split, instance_id_a, instance_id_b) VALUES (?, ?, ?, ?, ?)")
    BATCH = 100_000

    log("scanning DDIs for non-protected sources")
    membership, dropped, stats = [], [], Counter()
    inserted = 0
    for ddi_id, pfam_a, pfam_b, id_a, id_b, negative, source in external_ddis(conn):
        stats["ddis"] += 1
        if stats["ddis"] % 100_000 == 0:
            log(f"{stats['ddis']} external DDIs seen, {inserted + len(membership)} membership rows")
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
        if len(membership) >= BATCH:
            writer_cur.executemany(INSERT, membership)
            inserted += len(membership)
            membership.clear()

    writer_cur.executemany(INSERT, membership)
    inserted += len(membership)
    membership.clear()
    log(f"committing {inserted} membership rows")
    conn.commit()
    log("committed")

    with open(args.dropped_out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(DROPPED_COLUMNS)
        writer.writerows(sorted(dropped))

    for key in ("ddis", "pairs", "partial", "dropped_no_instances", "dropped_no_distinct_parent"):
        log(f"{key} = {stats[key]}")
    log(f"methods = {','.join(args.method)} split = {args.split}")
    conn.close()

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
