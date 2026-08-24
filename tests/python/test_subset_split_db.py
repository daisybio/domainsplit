#!/usr/bin/env python3
"""Checks for bin/subset_split_db.py.

``ddi_split_membership`` already says which DDIs belong to which split, so
subsetting is a pure SQL filter -- and the thing that can go wrong is the
cascade: keeping a split's DDIs but leaving another split's domains, instances,
proteins, GO terms or PPI edges behind would publish a database that leaks the
rest of the run.

Asserted here: only this split's DDIs survive; the domains, instances and
proteins reachable only from other splits are gone along with their GO and PPI
rows; and the membership table is narrowed to the one ``(method, split)`` so each
published DB describes only itself.

Run directly (`python3 tests/python/test_subset_split_db.py`) or via pytest.
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

SUBSET = os.path.join(BIN, "subset_split_db.py")

from ddi_db_utils import ensure_domains, insert_ddis  # noqa: E402
from domainsplit_schema import add_instance, ddi_rows, make_db  # noqa: E402

ENV = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))

SPLITS = {
    "train": ("PF00001", "PF00002"),
    "validation": ("PF00003", "PF00004"),
}


def build_db(path):
    make_db(path)
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA foreign_keys=ON")

    families = sorted({p for pair in SPLITS.values() for p in pair})
    ensure_domains(conn, families)
    insert_ddis(conn, list(SPLITS.values()), negative=False, source="3did")

    instance = {}
    for i, pfam in enumerate(families, start=1):
        instance[pfam] = add_instance(conn, pfam, f"P{i}", 10, 60)
        conn.execute("INSERT INTO domain_go_terms(domain_id, go_accession) "
                     "SELECT id, 'GO:0000001' FROM domain WHERE pfam_id = ?", (pfam,))
        conn.execute("INSERT INTO protein_go_terms(protein_id, go_accession) "
                     "SELECT id, 'GO:0000002' FROM protein WHERE uniprot_id = ?", (f"P{i}",))

    # One PPI edge inside the train split, one crossing into validation.
    conn.execute(
        "INSERT INTO protein_protein_interaction(protein_id_a, protein_id_b, score) "
        "SELECT a.id, b.id, 900 FROM protein a, protein b WHERE a.uniprot_id='P1' AND b.uniprot_id='P2'"
    )
    conn.execute(
        "INSERT INTO protein_protein_interaction(protein_id_a, protein_id_b, score) "
        "SELECT a.id, b.id, 800 FROM protein a, protein b WHERE a.uniprot_id='P1' AND b.uniprot_id='P3'"
    )

    for split, (pfam_a, pfam_b) in SPLITS.items():
        ddi_id = conn.execute(
            "SELECT ddi.id FROM domain_domain_interaction AS ddi "
            "JOIN domain AS da ON da.id = ddi.domain_id_a "
            "JOIN domain AS db ON db.id = ddi.domain_id_b "
            "WHERE da.pfam_id = ? AND db.pfam_id = ?", (pfam_a, pfam_b)
        ).fetchone()[0]
        conn.execute(
            "INSERT INTO ddi_split_membership(ddi_id, method, split, instance_id_a, instance_id_b) "
            "VALUES (?, 'minimal_leakage', ?, ?, ?)",
            (ddi_id, split, instance[pfam_a], instance[pfam_b]),
        )
    conn.commit()
    conn.close()


def test_subset_to_one_split():
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "minimal_leakage_train.sqlite3")
        build_db(db)

        subprocess.run(
            [sys.executable, SUBSET, "--db", db, "--method", "minimal_leakage", "--split", "train",
             "--versions", os.path.join(tmp, "versions.yml"), "--process-name", "TEST"],
            check=True, env=ENV, capture_output=True, text=True,
        )

        conn = sqlite3.connect(db)
        rows = ddi_rows(conn)
        domains = {r[0] for r in conn.execute("SELECT pfam_id FROM domain")}
        proteins = {r[0] for r in conn.execute("SELECT uniprot_id FROM protein")}
        instances = {r[0] for r in conn.execute("SELECT instance_id FROM domain_protein_map")}
        n_domain_go = conn.execute("SELECT COUNT(*) FROM domain_go_terms").fetchone()[0]
        n_protein_go = conn.execute("SELECT COUNT(*) FROM protein_go_terms").fetchone()[0]
        n_ppi = conn.execute("SELECT COUNT(*) FROM protein_protein_interaction").fetchone()[0]
        membership = conn.execute(
            "SELECT DISTINCT method, split FROM ddi_split_membership"
        ).fetchall()
        conn.close()

        assert rows == {("PF00001", "PF00002"): (0, "3did")}, rows
        assert domains == {"PF00001", "PF00002"}, domains
        assert proteins == {"P1", "P2"}, proteins
        assert instances == {"PF00001_P1_10_60", "PF00002_P2_10_60"}, instances
        assert (n_domain_go, n_protein_go) == (2, 2), (n_domain_go, n_protein_go)
        # The P1-P3 edge crossed into the validation split and went with P3.
        assert n_ppi == 1, n_ppi
        assert membership == [("minimal_leakage", "train")], membership


if __name__ == "__main__":
    test_subset_to_one_split()
    print("OK: subsetting keeps one split and cascades away everything else")
