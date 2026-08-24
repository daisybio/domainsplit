#!/usr/bin/env python3
"""Checks for bin/build_external_test.py.

Three properties:

  * only **external** DDIs are turned into a test set -- a 3did or
    sampled_negative pair is part of the splits, not of the held-out set;
  * every emitted instance pair has **distinct parent proteins**: two domains on
    one protein are not evidence that the two families interact, so a DDI whose
    only instances share a parent is dropped and reported;
  * sampling is deterministic under a fixed seed, and the seed is derived
    per-DDI, so the examples chosen for one DDI do not move when another DDI is
    added or removed.

Run directly (`python3 tests/python/test_build_external_test.py`) or via pytest.
"""

import csv
import os
import sqlite3
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

BUILDER = os.path.join(BIN, "build_external_test.py")

from ddi_db_utils import ensure_domains, insert_ddis, merge_external_ddis  # noqa: E402
from domainsplit_schema import add_instance, make_db  # noqa: E402

ENV = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))

METHODS = ["external_test", "external_test_hcni"]


def build_db(path):
    """One 3did pair, one rich external pair, and one external pair with a shared parent."""
    make_db(path)
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA foreign_keys=ON")

    ensure_domains(conn, {"PF00001", "PF00002"})
    insert_ddis(conn, [("PF00001", "PF00002")], negative=False, source="3did")
    merge_external_ddis(conn, [
        ("PF00010", "PF00011", 0, "PPIDM"),
        ("PF00020", "PF00021", 1, "negatome"),
    ])

    # UniProt accessions carry no underscore, so an instance id splits cleanly
    # into {family}_{uniprot}_{start}_{end}; keep the fixtures that way.
    add_instance(conn, "PF00001", "X1", 10, 60)
    add_instance(conn, "PF00002", "X2", 10, 60)
    # PF00010 x PF00011: three instances each, all on distinct proteins, so the
    # sampler has 9 candidate pairs to choose 3 from.
    for pfam, prefix in (("PF00010", "A"), ("PF00011", "B")):
        for i in (1, 2, 3):
            add_instance(conn, pfam, f"{prefix}{i}", 10 * i, 10 * i + 50)
    # PF00020 x PF00021: both families only occur on the same protein.
    add_instance(conn, "PF00020", "SHARED", 10, 60)
    add_instance(conn, "PF00021", "SHARED", 90, 140)
    conn.commit()
    conn.close()


def run_builder(db, tmp, seed=42, target=3, tag=""):
    dropped = os.path.join(tmp, f"dropped{tag}.tsv")
    args = [sys.executable, BUILDER, "--db", db]
    for method in METHODS:
        args += ["--method", method]
    args += ["--split", "test", "--target", str(target), "--pool-factor", "2",
             "--seed", str(seed), "--dropped-out", dropped,
             "--versions", os.path.join(tmp, f"versions{tag}.yml"), "--process-name", "TEST"]
    subprocess.run(args, check=True, env=ENV, capture_output=True, text=True)
    return dropped


def parents(db):
    """``{instance_id: uniprot_id}`` straight from domain_protein_map."""
    conn = sqlite3.connect(db)
    out = dict(conn.execute(
        "SELECT m.instance_id, p.uniprot_id FROM domain_protein_map AS m "
        "JOIN protein AS p ON p.id = m.protein_id"
    ))
    conn.close()
    return out


def membership(db):
    conn = sqlite3.connect(db)
    rows = conn.execute(
        "SELECT m.method, m.split, da.pfam_id, db.pfam_id, m.instance_id_a, m.instance_id_b "
        "FROM ddi_split_membership AS m "
        "JOIN domain_domain_interaction AS ddi ON ddi.id = m.ddi_id "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b "
        "ORDER BY m.method, m.instance_id_a, m.instance_id_b"
    ).fetchall()
    conn.close()
    return rows


def test_external_only_distinct_parents_and_report():
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "domainsplit.sqlite3")
        build_db(db)
        dropped_path = run_builder(db, tmp)

        rows = membership(db)
        families = {(r[2], r[3]) for r in rows}
        assert families == {("PF00010", "PF00011")}, families
        # Written under both external_test method directories, identically.
        by_method = {m: sorted((a, b) for meth, _, _, _, a, b in rows if meth == m) for m in METHODS}
        assert by_method[METHODS[0]] == by_method[METHODS[1]], by_method
        assert len(by_method[METHODS[0]]) == 3, by_method

        # Parents differ in every emitted pair, checked against the mapping table
        # rather than by parsing the instance id.
        parent = parents(db)
        for _, _, _, _, inst_a, inst_b in rows:
            assert parent[inst_a] != parent[inst_b], (inst_a, inst_b, parent[inst_a])

        with open(dropped_path) as fh:
            dropped = list(csv.DictReader(fh, delimiter="\t"))
        assert [(d["pfam_a"], d["pfam_b"], d["reason"]) for d in dropped] == [
            ("PF00020", "PF00021", "no_distinct_parent_pair")
        ], dropped


def test_seed_is_per_ddi_and_deterministic():
    with tempfile.TemporaryDirectory() as tmp:
        first = os.path.join(tmp, "a.sqlite3")
        build_db(first)
        run_builder(first, tmp, seed=42, tag="a")
        pairs_a = {(r[4], r[5]) for r in membership(first) if r[0] == METHODS[0]}

        # Same seed, same result.
        second = os.path.join(tmp, "b.sqlite3")
        build_db(second)
        run_builder(second, tmp, seed=42, tag="b")
        pairs_b = {(r[4], r[5]) for r in membership(second) if r[0] == METHODS[0]}
        assert pairs_a == pairs_b, (pairs_a, pairs_b)

        # Some other seed must move the sample -- otherwise --seed is not wired
        # in at all. Any single seed may coincide (3 pairs chosen from 9), so
        # the assertion is over a range.
        moved = False
        for seed in range(1, 12):
            other = os.path.join(tmp, f"c{seed}.sqlite3")
            build_db(other)
            run_builder(other, tmp, seed=seed, tag=f"c{seed}")
            if {(r[4], r[5]) for r in membership(other) if r[0] == METHODS[0]} != pairs_a:
                moved = True
                break
        assert moved, "no seed in 1..11 changed the sample; --seed is not reaching the sampler"


if __name__ == "__main__":
    test_external_only_distinct_parents_and_report()
    test_seed_is_per_ddi_and_deterministic()
    print("OK: external test set is external-only, parent-distinct and seed-deterministic")
