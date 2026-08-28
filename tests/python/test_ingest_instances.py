#!/usr/bin/env python3
"""Checks for bin/ingest_instances.py (no Nextflow, no cluster).

Three things are worth a test here:

  * ``start_pos``/``end_pos`` are stored as **integers**. This is the assertion
    whose absence let a real bug ship: the column was declared with no type, so
    it took BLOB (none) affinity and SQLite converted nothing on insert. This
    writer bound ``'10'`` while ENRICH's upsert bound ``10``, ``'10' != 10`` in a
    no-affinity column, so ``ON CONFLICT(domain_id, protein_id, start_pos,
    end_pos)`` never fired and every instance ended up as two rows -- the ingest
    row reachable by ``instance_id``, the enrich twin carrying the embeddings.
    Reading the columns back without checking their type could not see it;
  * the map is **instance level** -- a protein carrying two copies of one Pfam
    family must produce two ``domain_protein_map`` rows, because the embedding H5
    key contract (``{pfam}_{uniprot}_{start}_{end}``) and ppi-splitting's
    ``instances.tsv`` are both instance level. The old ``UNIQUE(domain_id,
    protein_id)`` silently kept one;
  * a non-human instance is a **hard failure**, not a filtered row: everything
    downstream (STRING, the UniProt idmapping) is human-only, so a
    non-9606 taxon means ppi-splitting was not run with
    ``instance_tier = 'human_only'`` and the resulting DB would be quietly wrong.

Run directly (`python3 tests/python/test_ingest_instances.py`) or via pytest.
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

INGEST = os.path.join(BIN, "ingest_instances.py")

from domainsplit_schema import make_db  # noqa: E402

INSTANCE_HEADER = ["instance_id", "family", "clan", "protein_id", "start", "end", "taxon_id", "source_db"]

# P11111 carries two copies of PF00001 -- the case the instance-level key exists for.
HUMAN_ROWS = [
    ("PF00001_P11111_10_60", "PF00001", "CL0001", "P11111", "10", "60", "9606", "reviewed"),
    ("PF00001_P11111_90_140", "PF00001", "CL0001", "P11111", "90", "140", "9606", "reviewed"),
    ("PF00002_Q22222_5_80", "PF00002", "CL0002", "Q22222", "5", "80", "9606", "unreviewed"),
]
MOUSE_ROW = ("PF00003_P33333_1_50", "PF00003", "CL0003", "P33333", "1", "50", "10090", "reviewed")


def write_instances(path, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(INSTANCE_HEADER) + "\n")
        for row in rows:
            fh.write("\t".join(row) + "\n")


def write_fasta(path, rows):
    with open(path, "w") as fh:
        for row in rows:
            fh.write(f">{row[0]}\nAAAA\n")


def run_ingest(tmp, rows, check=True):
    db = os.path.join(tmp, "domainsplit.sqlite3")
    make_db(db)
    instances = os.path.join(tmp, "instances.tsv")
    sequences = os.path.join(tmp, "sequences.fasta")
    mapping = os.path.join(tmp, "protein_domain_mapping.csv.gz")
    write_instances(instances, rows)
    write_fasta(sequences, rows)

    env = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))
    proc = subprocess.run(
        [sys.executable, INGEST, "--db", db, "--instances", instances,
         "--sequences", sequences, "--mapping-out", mapping,
         "--versions", os.path.join(tmp, "versions.yml"),
         "--process-name", "TEST:INGEST_INSTANCES"],
        check=check, env=env, capture_output=True, text=True,
    )
    return db, mapping, proc


def test_instance_level_rows():
    with tempfile.TemporaryDirectory() as tmp:
        db, mapping, _ = run_ingest(tmp, HUMAN_ROWS)

        conn = sqlite3.connect(db)
        rows = conn.execute(
            "SELECT d.pfam_id, p.uniprot_id, m.start_pos, m.end_pos, m.instance_id, m.clan, m.taxon_id "
            "FROM domain_protein_map AS m "
            "JOIN domain AS d ON d.id = m.domain_id "
            "JOIN protein AS p ON p.id = m.protein_id ORDER BY m.instance_id"
        ).fetchall()
        # Storage class, not just value: `typeof()` is what a no-affinity column
        # would answer 'text' to while `== 10` still looked right in Python.
        types = conn.execute(
            "SELECT DISTINCT typeof(start_pos), typeof(end_pos) FROM domain_protein_map"
        ).fetchall()
        n_domains = conn.execute("SELECT COUNT(*) FROM domain").fetchone()[0]
        n_proteins = conn.execute("SELECT COUNT(*) FROM protein").fetchone()[0]
        conn.close()

        assert len(rows) == 3, rows
        assert types == [("integer", "integer")], types
        for row in rows:
            assert isinstance(row[2], int), (row[2], type(row[2]))
            assert isinstance(row[3], int), (row[3], type(row[3]))
        # Both copies of PF00001 on P11111 survived as separate rows.
        assert [r[4] for r in rows] == [
            "PF00001_P11111_10_60", "PF00001_P11111_90_140", "PF00002_Q22222_5_80"
        ], rows
        assert n_domains == 2 and n_proteins == 2, (n_domains, n_proteins)

        with gzip.open(mapping, "rt") as fh:
            lines = fh.read().splitlines()
        assert lines[0] == "pfam_id,uniprot_id,start_pos,end_pos,sequence", lines[0]
        assert lines[1:] == [
            "PF00001,P11111,10,60,AAAA",
            "PF00001,P11111,90,140,AAAA",
            "PF00002,Q22222,5,80,AAAA",
        ], lines[1:]


def test_non_human_instance_is_fatal():
    with tempfile.TemporaryDirectory() as tmp:
        db, _, proc = run_ingest(tmp, HUMAN_ROWS + [MOUSE_ROW], check=False)

        assert proc.returncode != 0, proc.stdout
        assert "not human" in proc.stderr, proc.stderr
        assert "human_only" in proc.stderr, proc.stderr

        # Nothing was written -- the check runs before the first insert.
        conn = sqlite3.connect(db)
        assert conn.execute("SELECT COUNT(*) FROM domain_protein_map").fetchone()[0] == 0
        conn.close()


if __name__ == "__main__":
    test_instance_level_rows()
    test_non_human_instance_is_fatal()
    print("OK: ingest_instances is instance level and refuses non-human taxa")
