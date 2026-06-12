#!/usr/bin/env python3
"""Local integration check for the dual-source negative insertion (no cluster).

Validates the schema change (UNIQUE(domain_id_a, domain_id_b, source)) and
bin/insert_ppi_negative_selection.py together:

  * the four method labels are inserted, with 3did_random_addition copying the
    full 3did set and 3did_deletion only the pool-domain subset;
  * a pair can coexist under '3did' and '3did_random_addition' (duplicate by
    source);
  * the canonical sources still dedup across each other (a PPIDM pair equal to a
    3did pair is dropped) because insert_ddis defaults to dedup_across_sources.

Run directly or via pytest.
"""

import json
import os
import sqlite3
import subprocess
import sys
import tempfile

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)

from ddi_db_utils import ensure_domains, insert_ddis  # noqa: E402

INSERTER = os.path.join(BIN, "insert_ppi_negative_selection.py")

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


def pf(i):
    return f"PF{i:05d}"


def count(conn, source):
    return conn.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction WHERE source = ?",
        (source,),
    ).fetchone()[0]


def write_score(path, method):
    with open(path, "w") as fh:
        json.dump({
            "method": method, "seed": 7, "J": 0.1, "pa": 0.1, "deg": 0.1,
            "cov": 0.0, "n_sel": 2, "n_dom": 3, "mean_pa": 1.0,
            "pos_n_sel": 3, "pos_n_dom": 4, "pos_mean_pa": 2.0,
        }, fh)


def write_pairs(path, pairs):
    with open(path, "w") as fh:
        for a, b in pairs:
            fh.write(f"{a}\t{b}\n")


def test_dual_source_insert():
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "domainsplit.sqlite3")
        conn = sqlite3.connect(db)
        conn.executescript(SCHEMA)
        ensure_domains(conn, [pf(i) for i in range(1, 7)])

        # 3did positives, then a PPIDM batch overlapping (1,2).
        insert_ddis(conn, [(pf(1), pf(2)), (pf(1), pf(3)), (pf(1), pf(4))],
                    negative=False, source="3did")
        insert_ddis(conn, [(pf(1), pf(2)), (pf(5), pf(6))],
                    negative=False, source="PPIDM_Gold")
        conn.commit()
        assert count(conn, "3did") == 3
        assert count(conn, "PPIDM_Gold") == 1, "cross-source dedup broken"
        conn.close()

        # Pool covers only domains 1,2,3 -> 3did_deletion keeps (1,2),(1,3).
        pool = os.path.join(tmp, "neg_pool.npz")
        np.savez(pool, pool_dom=np.array([pf(1), pf(2), pf(3)], dtype=object))

        write_pairs(os.path.join(tmp, "pairs_deletion.tsv"),
                    [(pf(2), pf(5)), (pf(3), pf(6))])
        write_pairs(os.path.join(tmp, "pairs_random_addition.tsv"),
                    [(pf(1), pf(6)), (pf(4), pf(5))])
        write_score(os.path.join(tmp, "score_deletion.json"), "deletion")
        write_score(os.path.join(tmp, "score_random_addition.json"), "random_addition")

        env = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))
        subprocess.run(
            [sys.executable, INSERTER, "--db", db, "--pool", pool,
             "--pairs-deletion", os.path.join(tmp, "pairs_deletion.tsv"),
             "--pairs-random-addition", os.path.join(tmp, "pairs_random_addition.tsv"),
             "--score-deletion", os.path.join(tmp, "score_deletion.json"),
             "--score-random-addition", os.path.join(tmp, "score_random_addition.json"),
             "--scores-out", os.path.join(tmp, "scores.tsv")],
            check=True, env=env,
        )

        conn = sqlite3.connect(db)
        assert count(conn, "3did_random_addition") == 3, "full 3did copy wrong"
        assert count(conn, "3did_deletion") == 2, "pool-restricted copy wrong"
        assert count(conn, "inferred_ppi_screen_negative_for_deletion") == 2
        assert count(conn, "inferred_ppi_screen_negative_for_random_addition") == 2
        # original sources untouched
        assert count(conn, "3did") == 3
        assert count(conn, "PPIDM_Gold") == 1

        # The same pair (1,2) coexists under '3did' and '3did_random_addition'.
        n_dup = conn.execute(
            "SELECT COUNT(DISTINCT source) FROM domain_domain_interaction ddi "
            "JOIN domain da ON da.id = ddi.domain_id_a "
            "JOIN domain db ON db.id = ddi.domain_id_b "
            "WHERE da.pfam_id = ? AND db.pfam_id = ?",
            (pf(1), pf(2)),
        ).fetchone()[0]
        assert n_dup >= 2, f"(1,2) should exist under >=2 sources, got {n_dup}"
        conn.close()

    print("OK: dual-source insert + schema invariants hold")


if __name__ == "__main__":
    test_dual_source_insert()
