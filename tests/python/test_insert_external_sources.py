#!/usr/bin/env python3
"""End-to-end check for bin/insert_external_sources.py (no Nextflow, no cluster).

Seeds a database with a 3did pair, feeds the three normalized external TSVs the
PARSE_* modules produce, and asserts the whole contract at once:

  * the 3did pair keeps its own source and label, whatever the external sources
    claim about it, and that drop leaves no conflict-report noise;
  * agreeing external sources land on one row with a comma-joined source list;
  * disagreeing external sources take the pair out entirely and both appear in
    the conflict report.

Run directly (`python3 tests/python/test_insert_external_sources.py`) or via pytest.
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

INSERTER = os.path.join(BIN, "insert_external_sources.py")

from ddi_db_utils import ACTION_DROPPED, ensure_domains, insert_ddis  # noqa: E402
from domainsplit_schema import ddi_rows, make_db  # noqa: E402

HEADER = ["pfam_a", "pfam_b", "negative", "source"]

SINGLE_DOMAIN = [
    ("PF00001", "PF00002", 0, "single_domain_ppi"),  # collides with 3did
    ("PF00010", "PF00011", 0, "single_domain_ppi"),  # merges with PPIDM
]
PPIDM = [
    ("PF00010", "PF00011", 0, "PPIDM"),
    ("PF00010", "PF00011", 0, "PPIDM_Gold"),
    ("PF00020", "PF00021", 0, "PPIDM"),              # conflicts with negatome
    ("PF00020", "PF00021", 0, "PPIDM_Bronze"),
]
NEGATOME = [
    ("PF00020", "PF00021", 1, "negatome"),           # conflicts with PPIDM
    ("PF00030", "PF00031", 1, "negatome"),           # clean negative
]


def write_tsv(path, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(HEADER)
        writer.writerows(rows)


def test_insert_external_sources():
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "domainsplit.sqlite3")
        make_db(db)
        conn = sqlite3.connect(db)
        ensure_domains(conn, ["PF00001", "PF00002"])
        insert_ddis(conn, [("PF00001", "PF00002")], negative=False, source="3did")
        conn.commit()
        conn.close()

        paths = []
        for name, rows in (("single_domain_ppi", SINGLE_DOMAIN),
                           ("ppidm", PPIDM),
                           ("negatome", NEGATOME)):
            path = os.path.join(tmp, f"{name}_ddis.tsv")
            write_tsv(path, rows)
            paths.append(path)

        conflicts_out = os.path.join(tmp, "source_conflicts.tsv")
        env = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))
        subprocess.run(
            [sys.executable, INSERTER, "--db", db, "--ddis", *paths,
             "--conflicts-out", conflicts_out,
             "--versions", os.path.join(tmp, "versions.yml"),
             "--process-name", "TEST:INSERT_EXTERNAL_SOURCES"],
            check=True, env=env,
        )

        conn = sqlite3.connect(db)
        rows = ddi_rows(conn)
        conn.close()

        assert rows == {
            ("PF00001", "PF00002"): (0, "3did"),
            ("PF00010", "PF00011"): (0, "single_domain_ppi,PPIDM,PPIDM_Gold"),
            ("PF00030", "PF00031"): (1, "negatome"),
        }, rows

        with open(conflicts_out) as fh:
            conflicts = list(csv.DictReader(fh, delimiter="\t"))

        # Only the PPIDM/negatome disagreement is reported; the 3did collision is
        # ordinary and silent.
        assert {(c["pfam_a"], c["pfam_b"]) for c in conflicts} == {("PF00020", "PF00021")}
        assert sorted(c["source"] for c in conflicts) == [
            "PPIDM", "PPIDM_Bronze", "negatome"
        ]
        assert {c["action"] for c in conflicts} == {ACTION_DROPPED}


if __name__ == "__main__":
    test_insert_external_sources()
    print("OK: insert_external_sources merge/drop contract holds")
