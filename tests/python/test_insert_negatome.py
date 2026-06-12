#!/usr/bin/env python3
"""Local unit-check for bin/insert_negatome.py (no Nextflow, no cluster).

Builds a tiny empty Domainsplit SQLite and runs the Negatome inserter against a
small synthetic ``combined_pfam.txt``, asserting:

  * each whitespace-separated Pfam pair is stored as ``negative=1,
    source='negatome'`` with its domains auto-created;
  * lines without at least two tokens (blank / single-token) are skipped.

Run directly (`python3 tests/python/test_insert_negatome.py`) or via pytest.
"""

import os
import sqlite3
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
INSERTER = os.path.join(BIN, "insert_negatome.py")

SCHEMA = """
CREATE TABLE domain (id INTEGER PRIMARY KEY, pfam_id, name, UNIQUE(pfam_id));
CREATE TABLE domain_domain_interaction (
    id INTEGER PRIMARY KEY,
    domain_id_a, domain_id_b, negative,
    source VARCHAR(255),
    FOREIGN KEY(domain_id_a) REFERENCES domain ON DELETE CASCADE,
    FOREIGN KEY(domain_id_b) REFERENCES domain ON DELETE CASCADE,
    UNIQUE(domain_id_a, domain_id_b, source)
);
"""

NEGATOME_LINES = [
    "PF00001 PF00002",   # kept
    "PF00003\tPF00004",  # kept (tab separated)
    "PF00005 PF00006",   # kept
    "PF00007",           # single token -> skipped
    "",                  # blank -> skipped
]


def test_insert_negatome():
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "domainsplit.sqlite3")
        conn = sqlite3.connect(db)
        conn.executescript(SCHEMA)
        conn.commit()
        conn.close()

        negatome = os.path.join(tmp, "combined_pfam.txt")
        with open(negatome, "w") as fh:
            fh.write("\n".join(NEGATOME_LINES) + "\n")

        env = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))
        subprocess.run(
            [sys.executable, INSERTER, "--db", db, "--negatome", negatome,
             "--versions", os.path.join(tmp, "versions.yml"),
             "--process-name", "TEST:INSERT_NEGATOME"],
            check=True, env=env,
        )

        conn = sqlite3.connect(db)
        total = conn.execute("SELECT COUNT(*) FROM domain_domain_interaction").fetchone()[0]
        assert total == 3, f"expected 3 negatome DDIs, got {total}"

        rows = conn.execute(
            "SELECT COUNT(*) FROM domain_domain_interaction "
            "WHERE source = 'negatome' AND negative != 0"
        ).fetchone()[0]
        assert rows == 3, "all negatome rows must be negative with source 'negatome'"

        n_domains = conn.execute("SELECT COUNT(*) FROM domain").fetchone()[0]
        assert n_domains == 6, f"expected 6 auto-created domains, got {n_domains}"
        conn.close()


if __name__ == "__main__":
    test_insert_negatome()
    print("OK: insert_negatome invariants hold")
