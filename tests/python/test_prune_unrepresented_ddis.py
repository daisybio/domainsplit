#!/usr/bin/env python3
"""Checks for bin/prune_unrepresented_ddis.py.

The invariant this module owns is "every DDI in the master DB is representable
as an instance pair". A family can have zero instances -- a narrow ``--instance_tier``
drops families whose strata are all non-human, and a dead accession never
resolves -- and a DDI touching one of those cannot be split, embedded, or tested.

Asserted here: the DDI goes, whatever its source (3did included); the stranded
``domain`` and ``protein`` rows go with it; the split-membership rows follow
through ``ON DELETE CASCADE``; and every removal is named in the report.

Run directly (`python3 tests/python/test_prune_unrepresented_ddis.py`) or via pytest.
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

PRUNE = os.path.join(BIN, "prune_unrepresented_ddis.py")

from ddi_db_utils import ensure_domains, insert_ddis, merge_external_ddis  # noqa: E402
from domainsplit_schema import add_instance, ddi_rows, make_db  # noqa: E402

ENV = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))


def build_db(path):
    """PF00001/PF00002 are represented; PF00009 has no instance at all."""
    make_db(path)
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA foreign_keys=ON")

    ensure_domains(conn, {"PF00001", "PF00002", "PF00009"})
    insert_ddis(conn, [("PF00001", "PF00002"), ("PF00001", "PF00009")], negative=False, source="3did")
    merge_external_ddis(conn, [("PF00002", "PF00009", 0, "PPIDM")])

    inst_a = add_instance(conn, "PF00001", "P11111", 10, 60)
    inst_b = add_instance(conn, "PF00002", "Q22222", 5, 80)
    # An orphan protein with no domain instance of its own must go too.
    conn.execute("INSERT INTO protein(uniprot_id) VALUES ('P99999')")

    ddi_id = conn.execute(
        "SELECT ddi.id FROM domain_domain_interaction AS ddi "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b "
        "WHERE da.pfam_id = 'PF00001' AND db.pfam_id = 'PF00009'"
    ).fetchone()[0]
    conn.execute(
        "INSERT INTO ddi_split_membership(ddi_id, method, split, instance_id_a, instance_id_b) "
        "VALUES (?, 'minimal_leakage', 'train', ?, ?)",
        (ddi_id, inst_a, inst_b),
    )
    conn.commit()
    conn.close()


def test_prune():
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "domainsplit.sqlite3")
        report = os.path.join(tmp, "pruned_ddis.tsv")
        build_db(db)

        subprocess.run(
            [sys.executable, PRUNE, "--db", db, "--report-out", report,
             "--versions", os.path.join(tmp, "versions.yml"), "--process-name", "TEST"],
            check=True, env=ENV, capture_output=True, text=True,
        )

        conn = sqlite3.connect(db)
        rows = ddi_rows(conn)
        domains = {r[0] for r in conn.execute("SELECT pfam_id FROM domain")}
        proteins = {r[0] for r in conn.execute("SELECT uniprot_id FROM protein")}
        n_membership = conn.execute("SELECT COUNT(*) FROM ddi_split_membership").fetchone()[0]
        conn.close()

        assert rows == {("PF00001", "PF00002"): (0, "3did")}, rows
        assert domains == {"PF00001", "PF00002"}, domains
        assert proteins == {"P11111", "Q22222"}, proteins
        # The membership row pointed at a pruned DDI and cascaded away.
        assert n_membership == 0

        with open(report) as fh:
            pruned = list(csv.DictReader(fh, delimiter="\t"))
        assert {(p["pfam_a"], p["pfam_b"]) for p in pruned} == {
            ("PF00001", "PF00009"), ("PF00002", "PF00009")
        }, pruned
        # 3did is pruned like any other source -- being protected buys a pair
        # nothing once no instance can represent it.
        assert sorted(p["source"] for p in pruned) == ["3did", "PPIDM"], pruned
        assert all(p["reason"] == "no_instances:PF00009" for p in pruned), pruned


if __name__ == "__main__":
    test_prune()
    print("OK: unrepresentable DDIs, their domains, proteins and membership rows are pruned")
