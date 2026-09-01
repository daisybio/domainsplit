#!/usr/bin/env python3
"""Checks for bin/fetch_string_links.py (no Nextflow, no network).

The HTTP layer is stubbed with a fake ``requests.Session``, so what is under test
is the part that has to be right and cannot be seen from a finished run:

  * the taxon list comes from the **STRING id prefixes**, not from
    ``protein.taxon_id``. UniProt's ``OX`` can be a strain-level id with no STRING
    file behind it, so deriving from it would 404 on species STRING does cover;
  * a 404 is counted and reported, never fatal -- STRING does not cover every
    taxon UniProt does -- but every taxon 404ing *is* fatal, because a silently
    empty PPI table is what this whole path exists to prevent;
  * edges are filtered here, against the run's own proteins, so what reaches
    ``insert_ppi.py`` is small even when the universe is every species;
  * ``--min-proteins`` skips a taxon rather than paying a 50 MB download for one
    stray protein, and says so in the report.

Run directly (`python3 tests/python/test_fetch_string_links.py`) or via pytest.
"""

import gzip
import io
import os
import sqlite3
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pytest  # noqa: E402

import fetch_string_links as fsl  # noqa: E402
from domainsplit_schema import make_db  # noqa: E402

# Two human proteins, two mouse, one fly -- and the fly is the taxon STRING will
# not have a file for.
PROTEINS = [
    ("P11111", "9606", "9606.ENSP00000000001"),
    ("P22222", "9606", "9606.ENSP00000000002"),
    ("Q33333", "10090", "10090.ENSMUSP00000000001"),
    ("Q44444", "10090", "10090.ENSMUSP00000000002"),
    ("R55555", "7227", "7227.FBpp0000001"),
]

LINKS = {
    "9606": [
        ("9606.ENSP00000000001", "9606.ENSP00000000002", "900"),
        # Into a protein this run does not carry: dropped here, so insert_ppi
        # never sees it.
        ("9606.ENSP00000000001", "9606.ENSP00000999999", "700"),
    ],
    "10090": [("10090.ENSMUSP00000000001", "10090.ENSMUSP00000000002", "800")],
}


class FakeResponse:
    def __init__(self, status_code, payload=b""):
        self.status_code = status_code
        self._payload = payload

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def raise_for_status(self):
        if self.status_code >= 400:
            raise fsl.requests.HTTPError(f"status {self.status_code}")

    def iter_content(self, _chunk):
        yield self._payload


class FakeSession:
    """Serves LINKS; 404s for anything else. Records what was asked for."""

    def __init__(self, table=None):
        self.headers = {}
        self.requested = []
        self.table = LINKS if table is None else table

    def get(self, url, **_kwargs):
        self.requested.append(url)
        taxid = os.path.basename(url).split(".")[0]
        rows = self.table.get(taxid)
        if rows is None:
            return FakeResponse(404)
        buf = io.BytesIO()
        with gzip.GzipFile(fileobj=buf, mode="wb") as gz:
            gz.write(b"protein1 protein2 combined_score\n")
            for row in rows:
                gz.write((" ".join(row) + "\n").encode())
        return FakeResponse(200, buf.getvalue())


def build(tmp, proteins=PROTEINS):
    db = os.path.join(tmp, "domainsplit.sqlite3")
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
    return db, mapping


def run(tmp, argv_extra=(), session=None, monkeypatch=None, proteins=PROTEINS):
    db, mapping = build(tmp, proteins)
    out = os.path.join(tmp, "string_links.txt.gz")
    report = os.path.join(tmp, "string_taxa_report.tsv")
    session = session or FakeSession()
    monkeypatch.setattr(fsl.requests, "Session", lambda: session)
    monkeypatch.setattr(
        sys, "argv",
        ["fetch_string_links.py", "--db", db, "--string-map", mapping,
         "--out", out, "--report", report, "--jobs", "1",
         "--versions", os.path.join(tmp, "versions.yml"),
         "--process-name", "TEST:FETCH_STRING_LINKS", *argv_extra],
    )
    fsl.main()
    with gzip.open(out, "rt") as fh:
        links = [line.split() for line in fh.read().splitlines()]
    with open(report) as fh:
        rows = [line.split("\t") for line in fh.read().splitlines()]
    return session, links, {r[0]: dict(zip(rows[0], r)) for r in rows[1:]}


def test_taxa_come_from_the_string_id_prefix(monkeypatch, tmp_path):
    """`protein.taxon_id` is never read: the prefix is STRING's own species id."""
    session, links, report = run(str(tmp_path), monkeypatch=monkeypatch)

    asked = sorted(os.path.basename(u).split(".")[0] for u in session.requested)
    assert asked == ["10090", "7227", "9606"], asked
    # Kept edges only, both endpoints in the database, plus the header.
    assert links[0] == ["protein1", "protein2", "combined_score"]
    assert sorted(links[1:]) == sorted([
        ["9606.ENSP00000000001", "9606.ENSP00000000002", "900"],
        ["10090.ENSMUSP00000000001", "10090.ENSMUSP00000000002", "800"],
    ]), links
    assert report["9606"]["edges_read"] == "2"
    assert report["9606"]["edges_kept"] == "1"


def test_a_missing_taxon_is_reported_not_fatal():
    """STRING covers ~12k organisms and UniProt many more; a 404 is ordinary."""
    with tempfile.TemporaryDirectory() as tmp:
        with pytest.MonkeyPatch.context() as mp:
            _session, links, report = run(tmp, monkeypatch=mp)
    assert report["7227"]["status"] == "404"
    assert report["7227"]["edges_kept"] == "0"
    # The other two still contributed.
    assert len(links) == 3


def test_every_taxon_missing_is_fatal():
    """A green run with an empty PPI table is the failure this path exists to
    prevent, so 'nothing downloaded' has to stop the pipeline."""
    with tempfile.TemporaryDirectory() as tmp:
        with pytest.MonkeyPatch.context() as mp:
            with pytest.raises(SystemExit) as excinfo:
                run(tmp, session=FakeSession(table={}), monkeypatch=mp)
    assert "yielded a STRING links file" in str(excinfo.value)


def test_min_proteins_skips_a_thin_taxon():
    """One stray protein must not cost a 50 MB download -- but it has to be
    visible in the report, not silently absent."""
    with tempfile.TemporaryDirectory() as tmp:
        with pytest.MonkeyPatch.context() as mp:
            session, _links, report = run(
                tmp, argv_extra=["--min-proteins", "2"], monkeypatch=mp
            )
    asked = sorted(os.path.basename(u).split(".")[0] for u in session.requested)
    assert asked == ["10090", "9606"], asked
    assert report["7227"]["status"] == "skipped_min_proteins"


def test_no_string_xref_at_all_is_fatal():
    """An unreviewed-only universe carries almost no `DR STRING;` lines. That is a
    fact about UniProt, but a run with *none* would produce an empty PPI table, so
    it stops here with an explanation rather than there with nothing."""
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "domainsplit.sqlite3")
        make_db(db)
        conn = sqlite3.connect(db)
        conn.execute(
            "INSERT INTO protein(uniprot_id, taxon_id, reviewed) "
            "VALUES ('P11111', '9606', 'unreviewed')"
        )
        conn.commit()
        conn.close()
        mapping = os.path.join(tmp, "string_map.tsv.gz")
        with gzip.open(mapping, "wt") as fh:
            fh.write("")
        with pytest.MonkeyPatch.context() as mp:
            mp.setattr(fsl.requests, "Session", FakeSession)
            mp.setattr(
                sys, "argv",
                ["fetch_string_links.py", "--db", db, "--string-map", mapping,
                 "--out", os.path.join(tmp, "out.txt.gz"),
                 "--report", os.path.join(tmp, "report.tsv"),
                 "--versions", os.path.join(tmp, "versions.yml"),
                 "--process-name", "TEST"],
            )
            with pytest.raises(SystemExit) as excinfo:
                fsl.main()
    assert "carries a `DR STRING;` cross-reference" in str(excinfo.value)
