#!/usr/bin/env python3
"""INSERT_PROTEIN_SEQUENCES runs *after* ingest, so it must upsert, not insert.

Before the ppi-splitting integration, `protein` was empty when ENRICH ran, so the
script could insert blindly. Now INGEST_INSTANCES creates a bare `protein` row
per parent accession first, which changes what correct means: a plain INSERT
would hit UNIQUE(uniprot_id) on its first record and abort the run.

The other half of this file used to cover `insert_domain_protein_mapping.py`,
whose `ON CONFLICT` upsert existed only to attach ESM per-domain embeddings to
the rows ingest had already written. Embeddings left the database (they are
published as HDF5 by generate_domain_embeddings), which made that step a no-op
duplicate of ingest, so it is gone and so is its test. What replaced *its*
regression -- the upsert that silently never fired because
`domain_protein_map.start_pos` had no type affinity -- is asserted in
`test_ingest_instances.py`, at the writer, where a type can be checked.

Asserted here: the step adds no row, fills the column it owns, and leaves
`domain_protein_map` untouched.

Run directly (`python3 tests/python/test_enrich_after_ingest.py`) or via pytest.
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

INSERT_PROTEINS = os.path.join(BIN, "insert_protein_sequences.py")

from domainsplit_schema import add_instance, make_db  # noqa: E402

ENV = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))

# PF00001 sits twice on the same protein, which is the case the instance-level
# UNIQUE key exists for: a per-protein key would collapse these two rows into one.
INSTANCES = [
    ("PF00001", "P11111", 10, 20),
    ("PF00001", "P11111", 50, 60),
    ("PF00002", "Q22222", 5, 15),
]
SEQUENCES = {"P11111": "M" * 80, "Q22222": "K" * 40}


def build_inputs(tmp):
    """`(db, protein_domain_map, uniprot_fasta)` for the instances above."""
    db = os.path.join(tmp, "domainsplit.sqlite3")
    make_db(db)
    conn = sqlite3.connect(db)
    conn.execute("PRAGMA foreign_keys=ON")
    for pfam, uniprot, start, end in INSTANCES:
        add_instance(conn, pfam, uniprot, start, end, clan="CL0192", sequence="A" * (end - start + 1))
    conn.commit()
    conn.close()

    pd_map = os.path.join(tmp, "protein_domain_mapping.csv.gz")
    with gzip.open(pd_map, "wt", newline="") as fh:
        fh.write("pfam_id,uniprot_id,start_pos,end_pos,sequence\n")
        for pfam, uniprot, start, end in INSTANCES:
            fh.write(f"{pfam},{uniprot},{start},{end},{'A' * (end - start + 1)}\n")

    fasta = os.path.join(tmp, "uniprot_sequences.fasta.gz")
    with gzip.open(fasta, "wt") as fh:
        for uniprot, seq in SEQUENCES.items():
            fh.write(f">sp|{uniprot}|{uniprot}_HUMAN test protein\n{seq}\n")

    return db, pd_map, fasta


def run_insert_proteins(tmp, db, pd_map, fasta):
    """Run the real script the way the process does, and return the output DB."""
    # The script copies --db-in to ./domainsplit.sqlite3 in the cwd, so stage the
    # input under a different name.
    staged = os.path.join(tmp, "input.sqlite3")
    os.replace(db, staged)
    result = subprocess.run(
        [
            sys.executable, INSERT_PROTEINS,
            "--db-in", staged,
            "--uniprot-db", fasta,
            "--protein-domain-map", pd_map,
            "--versions", os.path.join(tmp, "v1.yml"),
            "--process-name", "TEST:INSERT_PROTEIN_SEQUENCES",
        ],
        env=ENV, cwd=tmp, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    return os.path.join(tmp, "domainsplit.sqlite3")


def test_protein_sequences_fill_ingested_rows_without_duplicating_them():
    with tempfile.TemporaryDirectory() as tmp:
        db, pd_map, fasta = build_inputs(tmp)
        out = run_insert_proteins(tmp, db, pd_map, fasta)

        conn = sqlite3.connect(out)
        assert conn.execute("SELECT COUNT(*) FROM protein").fetchone()[0] == len(SEQUENCES)
        # The step does not touch domain_protein_map at all any more.
        assert conn.execute("SELECT COUNT(*) FROM domain_protein_map").fetchone()[0] == len(INSTANCES)

        for uniprot, seq in SEQUENCES.items():
            row = conn.execute(
                "SELECT sequence FROM protein WHERE uniprot_id = ?", (uniprot,)
            ).fetchone()
            assert row[0] == seq, (uniprot, row)

        # The columns ingest owns survived, and the two copies of PF00001 are
        # still two rows.
        rows = conn.execute(
            "SELECT instance_id, clan, taxon_id FROM domain_protein_map ORDER BY start_pos"
        ).fetchall()
        assert len(rows) == len(INSTANCES)
        for instance_id, clan, taxon in rows:
            assert instance_id and instance_id.startswith("PF000")
            assert clan == "CL0192"
            assert str(taxon) == "9606"
        conn.close()


def test_no_embedding_columns_remain():
    """The schema mirror must not grow the columns back by accident."""
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "domainsplit.sqlite3")
        make_db(db)
        conn = sqlite3.connect(db)
        protein_cols = {r[1] for r in conn.execute("PRAGMA table_info(protein)")}
        map_cols = {r[1] for r in conn.execute("PRAGMA table_info(domain_protein_map)")}
        conn.close()
        assert not {c for c in protein_cols if "per_residue" in c}, protein_cols
        assert not {c for c in map_cols if "per_domain" in c}, map_cols


if __name__ == "__main__":
    test_protein_sequences_fill_ingested_rows_without_duplicating_them()
    test_no_embedding_columns_remain()
    print("ok")
