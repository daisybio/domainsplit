#!/usr/bin/env python3
"""Unit-check for the DDI merge/drop rules in bin/ddi_db_utils.py.

These rules are what make the source model work, so they are asserted directly
rather than through a module:

  * a pair held by a protected source (3did, sampled_negative) is never touched
    by an external source -- the incoming row is dropped silently, which is what
    guarantees the external test set is strictly unseen;
  * external sources that agree on ``negative`` merge into one row whose
    ``source`` is their comma-joined list;
  * external sources that disagree drop the pair entirely, deleting any existing
    external row for it, and report every dropped contribution;
  * pair order never matters -- (PF00002, PF00001) is the same pair as
    (PF00001, PF00002).

Run directly (`python3 tests/python/test_ddi_db_utils.py`) or via pytest.
"""

import os
import sqlite3
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ddi_db_utils import (  # noqa: E402
    ACTION_DELETED,
    ACTION_DROPPED,
    canonical_pair,
    count_source,
    ensure_domains,
    insert_ddis,
    merge_external_ddis,
)
from domainsplit_schema import ddi_rows, make_db  # noqa: E402


def _conn(tmp):
    db = os.path.join(tmp, "domainsplit.sqlite3")
    make_db(db)
    return sqlite3.connect(db)


def _seed_protected(conn, pairs, negative=False, source="3did"):
    ensure_domains(conn, [p for pair in pairs for p in pair])
    insert_ddis(conn, pairs, negative=negative, source=source)
    conn.commit()


def test_insert_into_empty():
    with tempfile.TemporaryDirectory() as tmp:
        conn = _conn(tmp)
        stats, conflicts = merge_external_ddis(
            conn, [("PF00002", "PF00001", 0, "single_domain_ppi")]
        )
        conn.commit()

        assert stats["inserted"] == 1
        assert conflicts == []
        # Stored in canonical (accession-number) order regardless of input order.
        assert ddi_rows(conn) == {("PF00001", "PF00002"): (0, "single_domain_ppi")}
        conn.close()


def test_protected_source_wins_silently():
    with tempfile.TemporaryDirectory() as tmp:
        conn = _conn(tmp)
        _seed_protected(conn, [("PF00001", "PF00002")])

        stats, conflicts = merge_external_ddis(
            conn,
            [
                ("PF00001", "PF00002", 0, "PPIDM"),      # same label
                ("PF00002", "PF00001", 1, "negatome"),   # opposite label, swapped
            ],
        )
        conn.commit()

        assert stats["dropped_protected"] == 1
        assert stats["conflict_pairs"] == 0, "a protected pair is not a conflict"
        assert conflicts == [], "protected drops are silent by design"
        assert ddi_rows(conn) == {("PF00001", "PF00002"): (0, "3did")}
        conn.close()


def test_external_sources_merge_on_agreement():
    with tempfile.TemporaryDirectory() as tmp:
        conn = _conn(tmp)
        stats, conflicts = merge_external_ddis(
            conn,
            [
                ("PF00001", "PF00002", 0, "single_domain_ppi"),
                ("PF00001", "PF00002", 0, "PPIDM"),
                ("PF00002", "PF00001", 0, "PPIDM_Gold"),  # swapped order, same pair
                ("PF00001", "PF00002", 0, "PPIDM"),       # exact duplicate
            ],
        )
        conn.commit()

        assert stats["inserted"] == 1
        assert conflicts == []
        rows = ddi_rows(conn)
        assert rows == {
            ("PF00001", "PF00002"): (0, "single_domain_ppi,PPIDM,PPIDM_Gold")
        }, rows
        conn.close()


def test_merge_appends_to_existing_external_row():
    """A later batch extends the source list instead of adding a second row."""
    with tempfile.TemporaryDirectory() as tmp:
        conn = _conn(tmp)
        merge_external_ddis(conn, [("PF00001", "PF00002", 0, "single_domain_ppi")])
        conn.commit()

        stats, conflicts = merge_external_ddis(
            conn,
            [
                ("PF00001", "PF00002", 0, "PPIDM"),
                ("PF00001", "PF00002", 0, "single_domain_ppi"),  # already there
            ],
        )
        conn.commit()

        assert stats["merged"] == 1 and stats["inserted"] == 0
        assert conflicts == []
        assert ddi_rows(conn) == {
            ("PF00001", "PF00002"): (0, "single_domain_ppi,PPIDM")
        }
        conn.close()


def test_label_conflict_within_one_batch_drops_the_pair():
    with tempfile.TemporaryDirectory() as tmp:
        conn = _conn(tmp)
        stats, conflicts = merge_external_ddis(
            conn,
            [
                ("PF00001", "PF00002", 0, "PPIDM"),
                ("PF00002", "PF00001", 1, "negatome"),
                ("PF00003", "PF00004", 0, "PPIDM"),  # untouched control
            ],
        )
        conn.commit()

        assert stats["conflict_pairs"] == 1
        assert stats["sources_dropped"] == 2
        assert ddi_rows(conn) == {("PF00003", "PF00004"): (0, "PPIDM")}
        assert sorted((c[2], c[3], c[4]) for c in conflicts) == [
            ("PPIDM", 0, ACTION_DROPPED),
            ("negatome", 1, ACTION_DROPPED),
        ]
        assert all(c[:2] == ("PF00001", "PF00002") for c in conflicts)
        conn.close()


def test_label_conflict_deletes_the_existing_external_row():
    with tempfile.TemporaryDirectory() as tmp:
        conn = _conn(tmp)
        merge_external_ddis(conn, [("PF00001", "PF00002", 0, "single_domain_ppi")])
        conn.commit()

        stats, conflicts = merge_external_ddis(
            conn, [("PF00001", "PF00002", 1, "negatome")]
        )
        conn.commit()

        assert stats["conflict_pairs"] == 1
        assert ddi_rows(conn) == {}, "the pair must be gone entirely"
        actions = sorted((c[2], c[4]) for c in conflicts)
        assert actions == [
            ("negatome", ACTION_DROPPED),
            ("single_domain_ppi", ACTION_DELETED),
        ], actions
        conn.close()


def test_conflict_is_not_resurrected_by_a_third_source():
    """Once a pair is poisoned in a batch, no other source in it revives the pair."""
    with tempfile.TemporaryDirectory() as tmp:
        conn = _conn(tmp)
        merge_external_ddis(
            conn,
            [
                ("PF00001", "PF00002", 0, "PPIDM"),
                ("PF00001", "PF00002", 1, "negatome"),
                ("PF00001", "PF00002", 0, "single_domain_ppi"),
            ],
        )
        conn.commit()
        assert ddi_rows(conn) == {}
        conn.close()


def test_count_source_matches_whole_tokens_only():
    with tempfile.TemporaryDirectory() as tmp:
        conn = _conn(tmp)
        merge_external_ddis(
            conn,
            [
                ("PF00001", "PF00002", 0, "PPIDM"),
                ("PF00001", "PF00002", 0, "PPIDM_Gold"),
            ],
        )
        conn.commit()

        assert count_source(conn, "PPIDM") == 1
        assert count_source(conn, "PPIDM_Gold") == 1
        assert count_source(conn, "PPIDM_Silver") == 0
        assert count_source(conn, "PPID") == 0, "must not match a source prefix"
        conn.close()


def test_canonical_pair_is_numeric():
    assert canonical_pair("PF00010", "PF00002") == ("PF00002", "PF00010")
    assert canonical_pair("PF00002", "PF00010") == ("PF00002", "PF00010")


if __name__ == "__main__":
    test_insert_into_empty()
    test_protected_source_wins_silently()
    test_external_sources_merge_on_agreement()
    test_merge_appends_to_existing_external_row()
    test_label_conflict_within_one_batch_drops_the_pair()
    test_label_conflict_deletes_the_existing_external_row()
    test_conflict_is_not_resurrected_by_a_third_source()
    test_count_source_matches_whole_tokens_only()
    test_canonical_pair_is_numeric()
    print("OK: ddi_db_utils merge/drop rules hold")
