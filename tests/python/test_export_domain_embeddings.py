#!/usr/bin/env python3
"""The published embedding HDF5 must satisfy the benchmark's own key contract.

This is the seam that fails *silently*. `daisybio-domainbenchmark` reads a domain
encoding as

    h5[pfam_id][COALESCE(instance_id, 'r' || rowid)]

and skips any instance pair it cannot resolve. A key layout that drifts therefore
does not raise: it yields zero training rows, which looks exactly like a model
that found no data. So the contract gets its own test rather than incidental
coverage, and the assertion is written as the benchmark's lookup, not as a
restatement of what the export happens to do.

No GPU, no model weights: the chunk is a fake HDF5 keyed by instance id, which is
all `run_embeddings.py` promises to produce.

Run directly (`python3 tests/python/test_export_domain_embeddings.py`) or via pytest.
"""

import os
import sqlite3
import subprocess
import sys
import tempfile

import h5py
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

EXPORT = os.path.join(BIN, "export_domain_embeddings.py")

from domainsplit_schema import add_instance, make_db  # noqa: E402

ENV = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))

DIM = 8
# PF00001 twice on one protein again: two instances of one family must land as two
# datasets under the *same* domain group, which a per-protein key would collapse.
INSTANCES = [
    ("PF00001", "P11111", 10, 20),
    ("PF00001", "P11111", 50, 60),
    ("PF00002", "Q22222", 5, 15),
    ("PF00003", "R33333", 1, 30),
]
# Dropped from the chunk on purpose: --embedding_max_len skips long sequences, so a
# missing vector is expected and must not fail the export.
UNEMBEDDED = "PF00003_R33333_1_30"


def build(tmp):
    db = os.path.join(tmp, "domainsplit.sqlite3")
    make_db(db)
    conn = sqlite3.connect(db)
    conn.execute("PRAGMA foreign_keys=ON")
    for pfam, uniprot, start, end in INSTANCES:
        add_instance(conn, pfam, uniprot, start, end)
    conn.commit()
    conn.close()

    # Two chunks, the way two FASTA shards arrive, staged as chunk1/chunk2.
    keys = [f"{p}_{u}_{s}_{e}" for p, u, s, e in INSTANCES if f"{p}_{u}_{s}_{e}" != UNEMBEDDED]
    for i, group in enumerate((keys[:2], keys[2:]), start=1):
        with h5py.File(os.path.join(tmp, f"chunk{i}"), "w") as h5:
            for j, key in enumerate(group):
                h5.create_dataset(key, data=np.full(DIM, float(j), dtype="float16"))
    return db


def run_export(tmp, db, model="esm3"):
    out = os.path.join(tmp, f"{model}_domain_embeddings.h5")
    done = subprocess.run(
        [sys.executable, EXPORT, "--db", db, "--chunk-glob", "chunk*",
         "--model", model, "--output-h5", out,
         "--versions", os.path.join(tmp, "versions.yml"), "--process-name", "TEST:EXPORT"],
        env=ENV, cwd=tmp, capture_output=True, text=True,
    )
    assert done.returncode == 0, done.stderr
    return out, done.stdout


def test_layout_matches_the_benchmark_lookup():
    with tempfile.TemporaryDirectory() as tmp:
        db = build(tmp)
        out, stdout = run_export(tmp, db)

        # The benchmark's own query, verbatim -- resolving domain_id to its Pfam
        # accession, which is what the outer HDF5 key is.
        conn = sqlite3.connect(db)
        rows = conn.execute(
            "SELECT d.pfam_id, COALESCE(m.instance_id, 'r' || m.rowid) "
            "FROM domain_protein_map AS m JOIN domain AS d ON d.id = m.domain_id"
        ).fetchall()
        conn.close()
        assert len(rows) == len(INSTANCES), rows

        with h5py.File(out, "r") as h5:
            resolved = []
            for pfam_id, key in rows:
                group = pfam_id
                if key == UNEMBEDDED:
                    # No vector was produced for it, so it must simply be absent --
                    # not present-but-empty, which would train on zeros.
                    assert group not in h5 or key not in h5[group], key
                    continue
                assert group in h5, f"no group {group!r} in {sorted(h5.keys())}"
                assert key in h5[group], f"{key!r} not in {sorted(h5[group].keys())}"
                assert h5[group][key].shape == (DIM,)
                assert h5[group][key].dtype == np.float16
                resolved.append(key)
        assert len(resolved) == len(INSTANCES) - 1, resolved

        # The two instances of PF00001 are two datasets in one group, and the
        # group is named by the accession -- not by a surrogate that changes run
        # to run, which is the whole point of the key.
        with h5py.File(out, "r") as h5:
            assert sorted(h5["PF00001"].keys()) == ["PF00001_P11111_10_20", "PF00001_P11111_50_60"]
            assert all(g.startswith("PF") for g in h5.keys()), sorted(h5.keys())

        assert "1 instances have no esm3 embedding" in stdout, stdout


def test_root_attributes_are_the_cross_run_guard():
    with tempfile.TemporaryDirectory() as tmp:
        db = build(tmp)
        out, _ = run_export(tmp, db, model="prott5")
        with h5py.File(out, "r") as h5:
            attrs = dict(h5.attrs)
        assert attrs["model"] == "prott5"
        assert attrs["pooling"] == "mean"
        assert attrs["dim"] == DIM
        assert attrs["dtype"] == "float16"
        assert attrs["key_layout"] == "{pfam_id}/{instance_id}"
        assert attrs["n_domains"] == 2, attrs          # PF00001 and PF00002
        assert attrs["n_instances"] == len(INSTANCES) - 1, attrs
        # No run identifier: with accession keys there is nothing for one to
        # guard, and a stale one would invite exactly the trust it cannot earn.
        assert "domainsplit_run" not in attrs, attrs


def test_total_key_mismatch_is_fatal():
    """An empty intersection is a wiring error, and must not publish a valid-looking file."""
    with tempfile.TemporaryDirectory() as tmp:
        db = build(tmp)
        with h5py.File(os.path.join(tmp, "chunk1"), "w") as h5:
            h5.create_dataset("PF99999_Z99999_1_2", data=np.zeros(DIM, dtype="float16"))
        os.remove(os.path.join(tmp, "chunk2"))
        done = subprocess.run(
            [sys.executable, EXPORT, "--db", db, "--chunk-glob", "chunk*", "--model", "esmc",
             "--output-h5", os.path.join(tmp, "esmc_domain_embeddings.h5"),
             "--versions", os.path.join(tmp, "v.yml"), "--process-name", "TEST"],
            env=ENV, cwd=tmp, capture_output=True, text=True,
        )
        assert done.returncode != 0, done.stdout
        assert "have diverged" in done.stderr + done.stdout


if __name__ == "__main__":
    test_layout_matches_the_benchmark_lookup()
    test_root_attributes_are_the_cross_run_guard()
    test_total_key_mismatch_is_fatal()
    print("OK: the exported layout satisfies the benchmark's key contract")
