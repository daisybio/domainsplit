#!/usr/bin/env python3
"""Local unit-check for bin/insert_ppidm.py (no Nextflow, no cluster).

Builds a tiny empty Domainsplit SQLite and runs the PPIDM inserter against a
small synthetic ``predicted_ddi_ppi.tsv``, asserting:

  * domain tokens like ``10114/PF00069`` are parsed down to the Pfam accession;
  * each kept row is stored as ``negative=0, source='PPIDM_<Class>'``;
  * classes are processed Gold -> Silver -> Bronze, so a pair appearing under
    two classes is kept only under the highest-confidence one (cross-source
    dedup in insert_ddis);
  * unparseable tokens are skipped, and ``--classes`` filters which classes are
    inserted at all.

Run directly (`python3 tests/python/test_insert_ppidm.py`) or via pytest.
"""

import os
import sqlite3
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
INSERTER = os.path.join(BIN, "insert_ppidm.py")

# Matches the schema the pipeline's INIT_DOMAINSPLIT_DB creates for these tables
# (see tests/python/test_insert_ppi_negative_selection.py).
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


def count(conn, source):
    return conn.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction WHERE source = ?",
        (source,),
    ).fetchone()[0]


def run_inserter(db, ppidm, classes, tmp):
    env = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))
    subprocess.run(
        [sys.executable, INSERTER, "--db", db, "--ppidm", ppidm,
         "--classes", classes,
         "--versions", os.path.join(tmp, "versions.yml"),
         "--process-name", "TEST:INSERT_PPIDM"],
        check=True, env=env,
    )


# Tokens carry a leading numeric id before the slash, as in real PPIDM output.
PPIDM_ROWS = [
    "domain_1\tdomain_2\tclass",          # header (skipped)
    "10/PF00001\t20/PF00002\tGold",       # kept -> PPIDM_Gold
    "30/PF00003\t40/PF00004\tSilver",     # kept -> PPIDM_Silver
    "50/PF00005\t60/PF00006\tBronze",     # kept -> PPIDM_Bronze
    "10/PF00001\t20/PF00002\tSilver",     # duplicate pair, lower class -> dropped
    "junk\tnonsense\tGold",               # unparseable -> skipped
]


def write_ppidm(path, rows):
    with open(path, "w") as fh:
        fh.write("\n".join(rows) + "\n")


def test_insert_ppidm_all_classes():
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "domainsplit.sqlite3")
        conn = sqlite3.connect(db)
        conn.executescript(SCHEMA)
        conn.commit()
        conn.close()

        ppidm = os.path.join(tmp, "predicted_ddi_ppi.tsv")
        write_ppidm(ppidm, PPIDM_ROWS)

        run_inserter(db, ppidm, "Bronze,Silver,Gold", tmp)

        conn = sqlite3.connect(db)
        # One pair per class; the duplicate (PF00001, PF00002) is kept only under
        # Gold (processed first) and dropped for Silver via cross-source dedup.
        assert count(conn, "PPIDM_Gold") == 1, "Gold count wrong"
        assert count(conn, "PPIDM_Silver") == 1, "Silver count wrong (dedup failed?)"
        assert count(conn, "PPIDM_Bronze") == 1, "Bronze count wrong"

        # All kept rows are positives stored under a PPIDM_* source only.
        total = conn.execute("SELECT COUNT(*) FROM domain_domain_interaction").fetchone()[0]
        assert total == 3, f"expected 3 DDIs total, got {total}"
        neg = conn.execute(
            "SELECT COUNT(*) FROM domain_domain_interaction WHERE negative != 0"
        ).fetchone()[0]
        assert neg == 0, "PPIDM rows must be positives"

        # The duplicate pair exists only under Gold, not Silver.
        n_sources = conn.execute(
            "SELECT COUNT(DISTINCT source) FROM domain_domain_interaction ddi "
            "JOIN domain da ON da.id = ddi.domain_id_a "
            "JOIN domain db ON db.id = ddi.domain_id_b "
            "WHERE da.pfam_id = ? AND db.pfam_id = ?",
            ("PF00001", "PF00002"),
        ).fetchone()[0]
        assert n_sources == 1, f"(PF00001,PF00002) should be under 1 source, got {n_sources}"
        conn.close()


def test_insert_ppidm_class_filter():
    """--classes restricts which classes are inserted at all."""
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "domainsplit.sqlite3")
        conn = sqlite3.connect(db)
        conn.executescript(SCHEMA)
        conn.commit()
        conn.close()

        ppidm = os.path.join(tmp, "predicted_ddi_ppi.tsv")
        write_ppidm(ppidm, PPIDM_ROWS)

        run_inserter(db, ppidm, "Gold", tmp)

        conn = sqlite3.connect(db)
        assert count(conn, "PPIDM_Gold") == 1
        assert count(conn, "PPIDM_Silver") == 0, "Silver should be excluded"
        assert count(conn, "PPIDM_Bronze") == 0, "Bronze should be excluded"
        conn.close()


if __name__ == "__main__":
    test_insert_ppidm_all_classes()
    test_insert_ppidm_class_filter()
    print("OK: insert_ppidm class handling + dedup invariants hold")
