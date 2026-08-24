#!/usr/bin/env python3
"""The two ENRICH steps that now run *after* ingest must upsert, not insert.

Before the ppi-splitting integration, `protein` and `domain_protein_map` were
empty when ENRICH ran, so both scripts could insert blindly. Now INGEST_INSTANCES
fills them from `instances.tsv` first, which changes what correct means:

- `insert_proteins_with_embeddings.py` would hit UNIQUE(uniprot_id) on its first
  record and abort the run;
- `insert_domain_protein_mapping.py` used `INSERT OR IGNORE`, so it would skip
  every existing row and the ESM per-domain embeddings would silently never be
  written -- the worse failure of the two, because nothing would report it.

Asserted here: neither step adds or duplicates a row, both fill the columns they
own, and the ingest-owned columns (instance_id, clan, taxon_id) survive.

Run directly (`python3 tests/python/test_enrich_after_ingest.py`) or via pytest.
"""

import gzip
import os
import pickle
import sqlite3
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

INSERT_PROTEINS = os.path.join(BIN, "insert_proteins_with_embeddings.py")
INSERT_MAPPING = os.path.join(BIN, "insert_domain_protein_mapping.py")

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
    import h5py
    import numpy as np

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

    prott5 = os.path.join(tmp, "prott5.h5")
    with h5py.File(prott5, "w") as h5:
        for uniprot, seq in SEQUENCES.items():
            h5.create_dataset(uniprot, data=np.zeros((len(seq), 4), dtype="float32"))

    esm_protein = os.path.join(tmp, "esm_protein.h5")
    with h5py.File(esm_protein, "w") as h5:
        for uniprot, seq in SEQUENCES.items():
            for model in ("esm3", "esmc"):
                h5.create_dataset(f"{uniprot}/{model}", data=np.ones((len(seq), 4), dtype="float32"))

    esm_domain = os.path.join(tmp, "esm_domain.h5")
    with h5py.File(esm_domain, "w") as h5:
        for pfam, uniprot, start, end in INSTANCES:
            for model in ("esm3", "esmc"):
                h5.create_dataset(
                    f"{pfam}_{uniprot}_{start}_{end}/{model}", data=np.full(4, 7.0, dtype="float32")
                )

    return db, pd_map, fasta, prott5, esm_protein, esm_domain


def run(cmd, tmp):
    result = subprocess.run(cmd, env=ENV, cwd=tmp, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return os.path.join(tmp, "domainsplit.sqlite3")


def test_enrich_fills_ingested_rows_without_duplicating_them():
    with tempfile.TemporaryDirectory() as tmp:
        db, pd_map, fasta, prott5, esm_protein, esm_domain = build_inputs(tmp)

        # Both scripts copy --db-in to ./domainsplit.sqlite3 in the cwd, so stage
        # the input under a different name and chain them the way ENRICH does.
        staged = os.path.join(tmp, "input.sqlite3")
        os.rename(db, staged)
        run(
            [
                sys.executable, INSERT_PROTEINS,
                "--db-in", staged,
                "--uniprot-db", fasta,
                "--protein-domain-map", pd_map,
                "--prott5-embeddings", prott5,
                "--esm-protein-embeddings", esm_protein,
                "--versions", os.path.join(tmp, "v1.yml"),
                "--process-name", "TEST:INSERT_PROTEINS",
            ],
            tmp,
        )
        os.replace(os.path.join(tmp, "domainsplit.sqlite3"), staged)
        out = run(
            [
                sys.executable, INSERT_MAPPING,
                "--db-in", staged,
                "--protein-domain-map", pd_map,
                "--esm-domain-embeddings", esm_domain,
                "--versions", os.path.join(tmp, "v2.yml"),
                "--process-name", "TEST:INSERT_MAPPING",
            ],
            tmp,
        )

        conn = sqlite3.connect(out)
        assert conn.execute("SELECT COUNT(*) FROM protein").fetchone()[0] == len(SEQUENCES)
        assert conn.execute("SELECT COUNT(*) FROM domain_protein_map").fetchone()[0] == len(INSTANCES)

        for uniprot, seq in SEQUENCES.items():
            row = conn.execute(
                "SELECT sequence, prott5_per_residue, esm3_per_residue, esmc_per_residue "
                "FROM protein WHERE uniprot_id = ?",
                (uniprot,),
            ).fetchone()
            assert row[0] == seq
            assert all(value is not None for value in row[1:])

        # Every instance got its per-domain embeddings, and the columns ingest
        # owns were not overwritten on the way.
        rows = conn.execute(
            "SELECT instance_id, clan, taxon_id, esm3_per_domain, esmc_per_domain "
            "FROM domain_protein_map ORDER BY start_pos"
        ).fetchall()
        assert len(rows) == len(INSTANCES)
        for instance_id, clan, taxon, esm3, esmc in rows:
            assert instance_id and instance_id.startswith("PF000")
            assert clan == "CL0192"
            assert str(taxon) == "9606"
            assert pickle.loads(esm3)[0] == 7.0
            assert pickle.loads(esmc)[0] == 7.0
        conn.close()


if __name__ == "__main__":
    test_enrich_fills_ingested_rows_without_duplicating_them()
    print("ok")
