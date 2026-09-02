#!/usr/bin/env python3
"""Checks for bin/insert_ppi.py (no Nextflow, no cluster).

This script used to be ``pd.read_csv(links_file)`` -- fine for one 70 MB human
file, impossible for FETCH_STRING_LINKS' concatenation of every taxon in the run.
Three things are worth asserting about the streaming rewrite:

  * it inserts the same edges the DataFrame version did, from a **multi-species**
    file, with no taxon bookkeeping -- STRING ids are globally unique because they
    carry their own taxon prefix;
  * the drop classes are counted **per taxon**. With one species in the file the
    totals were enough; with thousands, "none of these edges resolved" has to say
    *which* species, or a single misconfigured organism is invisible inside a
    large ``not_in_db``;
  * a file that resolves to nothing at all is still fatal. That guard is the whole
    reason the counters exist: before it, the wrong organism's links file produced
    zero PPI rows and a green run.

Run directly (`python3 tests/python/test_insert_ppi.py`) or via pytest.
"""

import gzip
import os
import sqlite3
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

INSERT_PPI = os.path.join(BIN, "insert_ppi.py")

from domainsplit_schema import make_db  # noqa: E402

ENV = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))

# Two species. The mouse pair is the one a human-only pipeline would have lost.
PROTEINS = [
    ("P11111", "9606", "9606.ENSP00000000001"),
    ("P22222", "9606", "9606.ENSP00000000002"),
    ("Q33333", "10090", "10090.ENSMUSP00000000001"),
    ("Q44444", "10090", "10090.ENSMUSP00000000002"),
]


def build(tmp, links_rows, proteins=PROTEINS):
    db = os.path.join(tmp, "in.sqlite3")
    make_db(db)
    conn = sqlite3.connect(db)
    conn.executemany(
        "INSERT INTO protein(uniprot_id, taxon_id, reviewed) VALUES (?, ?, 'reviewed')",
        [(acc, taxon) for acc, taxon, _sid in proteins],
    )
    conn.commit()
    conn.close()

    mapping = os.path.join(tmp, "string_map.tsv.gz")
    with gzip.open(mapping, "wt") as fh:
        for acc, _taxon, sid in proteins:
            fh.write(f"{acc}\tSTRING\t{sid}\n")

    links = os.path.join(tmp, "links.txt.gz")
    with gzip.open(links, "wt") as fh:
        fh.write("protein1 protein2 combined_score\n")
        for row in links_rows:
            fh.write(" ".join(row) + "\n")
    return db, mapping, links


def run(tmp, links_rows, proteins=PROTEINS, check=True):
    db, mapping, links = build(tmp, links_rows, proteins)
    proc = subprocess.run(
        [sys.executable, INSERT_PPI, "--db-in", db, "--string-ppi", links,
         "--uniprot-id-mapping", mapping,
         "--versions", os.path.join(tmp, "versions.yml"),
         "--process-name", "TEST:INSERT_PPI"],
        cwd=tmp, check=check, env=ENV, capture_output=True, text=True,
    )
    return proc


def edges(tmp):
    conn = sqlite3.connect(os.path.join(tmp, "domainsplit.sqlite3"))
    rows = set(
        conn.execute(
            "SELECT a.uniprot_id, b.uniprot_id, ppi.score "
            "FROM protein_protein_interaction AS ppi "
            "JOIN protein AS a ON a.id = ppi.protein_id_a "
            "JOIN protein AS b ON b.id = ppi.protein_id_b"
        )
    )
    conn.close()
    return rows


def test_multi_species_links_insert_without_taxon_bookkeeping():
    rows = [
        ("9606.ENSP00000000001", "9606.ENSP00000000002", "900"),
        ("10090.ENSMUSP00000000001", "10090.ENSMUSP00000000002", "800"),
        # An edge into a protein the run does not carry: ordinary, and the large
        # number on a real file.
        ("9606.ENSP00000000001", "9606.ENSP00000999999", "700"),
    ]
    with tempfile.TemporaryDirectory() as tmp:
        proc = run(tmp, rows)
        got = edges(tmp)

    assert got == {
        ("P11111", "P22222", 900.0),
        ("Q33333", "Q44444", 800.0),
    }, got
    assert "3 read, 2 inserted" in proc.stdout, proc.stdout


def test_per_taxon_counters_localise_a_broken_join():
    """One species resolving to nothing must be named, not buried in the total."""
    rows = [
        ("9606.ENSP00000000001", "9606.ENSP00000000002", "900"),
        # Mouse ids the mapping does not know: a wrong-release or wrong-organism
        # links file for that species alone.
        ("10090.ENSMUSP00000009998", "10090.ENSMUSP00000009999", "800"),
        ("10090.ENSMUSP00000009997", "10090.ENSMUSP00000009999", "700"),
    ]
    with tempfile.TemporaryDirectory() as tmp:
        proc = run(tmp, rows)

    assert "9606: read=1 inserted=1" in proc.stdout, proc.stdout
    assert "10090: read=2 inserted=0 unmapped=2" in proc.stdout, proc.stdout
    assert "1 taxa contributed edges but none that resolved" in proc.stdout, proc.stdout


def test_a_file_that_resolves_to_nothing_is_fatal():
    rows = [("7227.FBpp0000001", "7227.FBpp0000002", "900")]
    with tempfile.TemporaryDirectory() as tmp:
        proc = run(tmp, rows, check=False)

    assert proc.returncode != 0
    assert "none of 1 STRING edges could be inserted" in proc.stderr, proc.stderr


def test_duplicate_edges_do_not_abort_the_insert():
    """A concatenation of per-organism files can repeat a cross-species edge, and
    ``UNIQUE(protein_id_a, protein_id_b)`` would otherwise take the run down at the
    end of a long insert."""
    rows = [
        ("9606.ENSP00000000001", "9606.ENSP00000000002", "900"),
        ("9606.ENSP00000000001", "9606.ENSP00000000002", "900"),
    ]
    with tempfile.TemporaryDirectory() as tmp:
        run(tmp, rows)
        assert edges(tmp) == {("P11111", "P22222", 900.0)}


def test_a_non_numeric_score_is_fatal():
    """`score` reaches consumers as a number or the run stops here.

    It used to reach them as whatever string the links file held: the column is
    typeless in the schema, INSERT bound `parts[2]` verbatim, and a downstream
    numeric confidence filter (`score >= 400`) then raised a TypeError deep in a
    benchmark run instead of failing at ingestion. A third field that is not a
    number means this is not a STRING links file, so say so here.
    """
    rows = [("9606.ENSP00000000001", "9606.ENSP00000000002", "not-a-score")]
    with tempfile.TemporaryDirectory() as tmp:
        proc = run(tmp, rows, check=False)

    assert proc.returncode != 0
    assert "non-numeric combined_score" in proc.stderr, proc.stderr


def test_scores_are_stored_as_numbers():
    """The stored class, not just the value: REAL affinity plus a parsed insert."""
    rows = [("9606.ENSP00000000001", "9606.ENSP00000000002", "900")]
    with tempfile.TemporaryDirectory() as tmp:
        run(tmp, rows)
        conn = sqlite3.connect(os.path.join(tmp, "domainsplit.sqlite3"))
        classes = {
            r[0]
            for r in conn.execute(
                "SELECT typeof(score) FROM protein_protein_interaction"
            )
        }
        conn.close()

    assert classes == {"real"}, classes


if __name__ == "__main__":
    test_multi_species_links_insert_without_taxon_bookkeeping()
    test_per_taxon_counters_localise_a_broken_join()
    test_a_file_that_resolves_to_nothing_is_fatal()
    test_a_non_numeric_score_is_fatal()
    test_scores_are_stored_as_numbers()
    test_duplicate_edges_do_not_abort_the_insert()
    print("ok")
