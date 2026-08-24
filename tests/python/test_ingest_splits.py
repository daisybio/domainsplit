#!/usr/bin/env python3
"""Checks for bin/ingest_sampled_negatives.py and bin/ingest_split_membership.py.

Covers the three properties the ingest boundary has to hold:

  * a ``label=0`` family pair becomes exactly one ``sampled_negative`` DDI no
    matter how many splits or negative sets produced it;
  * ``val`` is renamed to ``validation`` here, and nowhere else;
  * running before INSERT_EXTERNAL_SOURCES is not a suggestion -- if an external
    source already holds a sampled-negative pair, the insertion order that makes
    the external test set unseen has been broken, and the script must say so
    rather than carry on.

Run directly (`python3 tests/python/test_ingest_splits.py`) or via pytest.
"""

import os
import sqlite3
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

NEGATIVES = os.path.join(BIN, "ingest_sampled_negatives.py")
MEMBERSHIP = os.path.join(BIN, "ingest_split_membership.py")

from ddi_db_utils import ensure_domains, insert_ddis, merge_external_ddis  # noqa: E402
from domainsplit_schema import add_instance, ddi_rows, make_db  # noqa: E402

ENV = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))

# Two 3did positives and, in the splits below, one sampled negative per split.
POSITIVES = [("PF00001", "PF00002"), ("PF00003", "PF00004")]

FAMILY_SPLITS = {
    "train": [("PF00001", "PF00002", 1), ("PF00001", "PF00003", 0)],
    "val": [("PF00003", "PF00004", 1), ("PF00002", "PF00004", 0)],
}
# The same negative pair reappears in the second negative set; it must not
# produce a second DDI row.
FAMILY_SPLITS_HCNI = {
    "train": [("PF00001", "PF00002", 1), ("PF00001", "PF00003", 0)],
}


def write_family_csv(path, rows):
    with open(path, "w") as fh:
        fh.write("protein1,protein2,label\n")
        for a, b, label in rows:
            fh.write(f"{a},{b},{label}\n")


def write_instance_csv(path, rows):
    with open(path, "w") as fh:
        fh.write("protein1,protein2,label,family1,family2\n")
        for ia, ib, label, fa, fb in rows:
            fh.write(f"{ia},{ib},{label},{fa},{fb}\n")


def seed_db(path):
    """3did positives plus one instance per family, as INGEST_INSTANCES would leave it."""
    make_db(path)
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA foreign_keys=ON")
    ensure_domains(conn, {p for pair in POSITIVES for p in pair})
    insert_ddis(conn, POSITIVES, negative=False, source="3did")
    instances = {}
    for i, pfam in enumerate(["PF00001", "PF00002", "PF00003", "PF00004"], start=1):
        instances[pfam] = add_instance(conn, pfam, f"P{i}{i}{i}{i}{i}{i}", 10, 60)
    conn.commit()
    conn.close()
    return instances


def run(script, db, specs, extra=(), check=True):
    args = [sys.executable, script, "--db", db]
    for spec in specs:
        args += ["--split", spec]
    args += list(extra)
    args += ["--versions", db + ".versions.yml", "--process-name", "TEST"]
    return subprocess.run(args, check=check, env=ENV, capture_output=True, text=True)


def test_sampled_negatives_dedup_and_rename():
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "domainsplit.sqlite3")
        instances = seed_db(db)

        specs = []
        for method, splits in (("minimal_leakage", FAMILY_SPLITS),
                               ("minimal_leakage_hcni", FAMILY_SPLITS_HCNI)):
            for split, rows in splits.items():
                path = os.path.join(tmp, f"{method}_{split}.csv")
                write_family_csv(path, rows)
                specs.append(f"{method}:{split}:{path}")
        run(NEGATIVES, db, specs)

        conn = sqlite3.connect(db)
        rows = ddi_rows(conn)
        conn.close()

        # PF00001/PF00003 appeared in two negative sets and is still one row.
        assert rows == {
            ("PF00001", "PF00002"): (0, "3did"),
            ("PF00003", "PF00004"): (0, "3did"),
            ("PF00001", "PF00003"): (1, "sampled_negative"),
            ("PF00002", "PF00004"): (1, "sampled_negative"),
        }, rows

        # Membership, over the same splits at instance level.
        specs = []
        for method, splits in (("minimal_leakage", FAMILY_SPLITS),
                               ("minimal_leakage_hcni", FAMILY_SPLITS_HCNI)):
            for split, fam_rows in splits.items():
                path = os.path.join(tmp, f"{method}_{split}_instances.csv")
                write_instance_csv(path, [
                    (instances[a], instances[b], label, a, b) for a, b, label in fam_rows
                ])
                specs.append(f"{method}:{split}:{path}")
        unresolved = os.path.join(tmp, "unresolved.tsv")
        run(MEMBERSHIP, db, specs, extra=["--unresolved-out", unresolved])

        conn = sqlite3.connect(db)
        seen = conn.execute(
            "SELECT method, split, COUNT(*) FROM ddi_split_membership "
            "GROUP BY method, split ORDER BY method, split"
        ).fetchall()
        conn.close()

        # `val` became `validation` at this boundary, and only here.
        assert seen == [
            ("minimal_leakage", "train", 2),
            ("minimal_leakage", "validation", 2),
            ("minimal_leakage_hcni", "train", 2),
        ], seen

        with open(unresolved) as fh:
            assert len(fh.read().splitlines()) == 1, "expected header only"


def test_external_source_before_negatives_is_fatal():
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "domainsplit.sqlite3")
        seed_db(db)

        # INSERT_EXTERNAL_SOURCES ran too early and claimed a pair the sampler
        # is about to produce.
        conn = sqlite3.connect(db)
        conn.execute("PRAGMA foreign_keys=ON")
        merge_external_ddis(conn, [("PF00001", "PF00003", 0, "PPIDM")])
        conn.commit()
        conn.close()

        path = os.path.join(tmp, "train.csv")
        write_family_csv(path, FAMILY_SPLITS["train"])
        proc = run(NEGATIVES, db, [f"minimal_leakage:train:{path}"], check=False)

        assert proc.returncode != 0, proc.stdout
        assert "already held by an external source" in proc.stderr, proc.stderr


if __name__ == "__main__":
    test_sampled_negatives_dedup_and_rename()
    test_external_source_before_negatives_is_fatal()
    print("OK: sampled negatives dedup, val->validation rename, insertion order enforced")
