#!/usr/bin/env python3
"""ENRICH_STRUCTURAL's AF confidence and missing-structure reports.

Asserted here:

  * the PAE submatrix stored at pae/<group>/<key> is the model's PAE over the
    two sliced windows in the *row's* orientation ([domain A, domain B]), for
    both an AFDB heterodimer (PAE rows in PDB chain/residue order) and an own
    AF3 run (rows mapped through token_chain_ids/token_res_ids), and AF3's
    chain_pair_iptm is reordered the same way;
  * a missing confidence file leaves NaN attrs / no PAE and is not fatal; a
    present PAE whose shape does not fit the model is fatal;
  * structure_coverage.tsv and af3_to_model/*.tsv list every membership row
    without structure with its reason.
"""

import csv
import json
import math
import os
import subprocess
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from test_structure_scoring import BIN, ENV, N_RES, build, pdb_text  # noqa: E402

# Value = 40*row + col: every entry is distinct and exact in float16.
PAE = (np.arange(2 * N_RES)[:, None] * 40 + np.arange(2 * N_RES)[None, :]).astype(float)
# Row protein P1 is file chain B (offset N_RES), window 5-15; P2 is chain A, window 3-9.
ROW_IDX = [N_RES + r - 1 for r in range(5, 16)] + [r - 1 for r in range(3, 10)]


def write_meta(tmp, name, source):
    with open(os.path.join(tmp, "af_metadata.csv"), "w") as fh:
        fh.write("uniprot_id_a,uniprot_id_b,path,source\n")
        fh.write(f"P2,P1,{name},{source}\n")


def afdb_fixture(tmp, pae=PAE, write_pae=True):
    import zstandard

    db, _, ia, ib, ic = build(tmp)
    os.rename(os.path.join(tmp, "p2_p1.pdb.zst"), os.path.join(tmp, "AF-P2-P1-model_v1.pdb.zst"))
    write_meta(tmp, "AF-P2-P1-model_v1.pdb.zst", "downloaded")
    if write_pae:
        payload = json.dumps([{"predicted_aligned_error": pae.tolist(),
                               "max_predicted_aligned_error": 31.75}]).encode()
        with open(os.path.join(tmp, "AF-P2-P1-predicted_aligned_error_v1.json.zst"), "wb") as fh:
            fh.write(zstandard.ZstdCompressor().compress(payload))
    return db, ia, ib, ic


def af3_fixture(tmp):
    import io

    from Bio.PDB.mmcifio import MMCIFIO
    from Bio.PDB.PDBParser import PDBParser

    db, _, ia, ib, ic = build(tmp)
    run = os.path.join(tmp, "p2_p1")
    os.makedirs(run)
    s = PDBParser(QUIET=True).get_structure("c", io.StringIO(pdb_text([("A", 4.0), ("B", 0.0)])))
    cif = MMCIFIO()
    cif.set_structure(s)
    cif.save(os.path.join(run, "p2_p1_model.cif"))
    with open(os.path.join(run, "p2_p1_summary_confidences.json"), "w") as fh:
        json.dump({"iptm": 0.81, "ptm": 0.77, "ranking_score": 0.9,
                   "chain_pair_iptm": [[0.9, 0.7], [0.6, 0.95]]}, fh)
    with open(os.path.join(run, "p2_p1_confidences.json"), "w") as fh:
        json.dump({"pae": PAE.tolist(),
                   "token_chain_ids": ["A"] * N_RES + ["B"] * N_RES,
                   "token_res_ids": list(range(1, N_RES + 1)) * 2}, fh)
    write_meta(tmp, "p2_p1/p2_p1_model.cif", "inferred")
    return db, ia, ib, ic


def enrich(tmp, db, check=True):
    return subprocess.run(
        [sys.executable, os.path.join(BIN, "enrich_structural_af.py"),
         "--db_in", db, "--db_out", os.path.join(tmp, "out.sqlite3"),
         "--structures_h5", os.path.join(tmp, "structures.h5"),
         "--af_metadata", os.path.join(tmp, "af_metadata.csv"), "--methods", "random,random_ilp",
         "--cpus", "1", "--max_mismatch_rate", "0.02",
         "--coverage", os.path.join(tmp, "coverage.tsv"), "--af3_to_model_dir", os.path.join(tmp, "af3"),
         "--versions", os.path.join(tmp, "versions.yml"), "--process_name", "TEST"],
        check=check, env=ENV, capture_output=True, text=True,
    )


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def stored(tmp, ia, ib):
    import h5py

    import utils_struct

    group, key = utils_struct.structure_h5_keys(ia, ib)
    with h5py.File(os.path.join(tmp, "structures.h5"), "r") as h5:
        attrs = dict(h5[group][key].attrs)
        pae = h5.get(f"pae/{group}/{key}")
        pae = None if pae is None else (pae[()], dict(pae.attrs))
        assert utils_struct.get_structure_bytes(h5, ia, ib) is not None
    return attrs, pae


def test_afdb_pae_is_sliced_in_row_orientation_and_reports_list_missing_pairs():
    with tempfile.TemporaryDirectory() as tmp:
        db, ia, ib, ic = afdb_fixture(tmp)
        log = enrich(tmp, db).stdout
        attrs, (pae, pattrs) = stored(tmp, ia, ib)
        assert attrs["af_source"] == "downloaded"
        assert all(math.isnan(attrs[k]) for k in ("iptm", "ptm", "ranking_score"))
        assert "chain_pair_iptm" not in attrs
        assert pae.dtype == np.float32
        np.testing.assert_array_equal(pae, PAE[np.ix_(ROW_IDX, ROW_IDX)])
        assert list(pattrs["residues_a"]) == list(range(5, 16))
        assert list(pattrs["residues_b"]) == list(range(3, 10))

        cov = {(r["method"], r["split"], r["label"]): r for r in read_tsv(os.path.join(tmp, "coverage.tsv"))}
        assert cov[("random", "train", "positive")]["with_pae"] == "1"
        assert cov[("random", "train", "negative")]["with_structure"] == "0"
        assert cov[("random", "train", "negative")]["fraction"] == "0.0"
        assert ("random_ilp", "train", "positive") in cov

        todo = read_tsv(os.path.join(tmp, "af3", "random.tsv"))
        assert todo == [{"uniprot_id_a": "P1", "uniprot_id_b": "P3", "reason": "no_af_metadata_entry",
                         "n_instance_pairs": "1", "labels": "negative", "splits": "train"}]
        assert read_tsv(os.path.join(tmp, "af3", "random_ilp.tsv")) == []
        assert read_tsv(os.path.join(tmp, "af3", "all_methods.tsv"))[0]["methods"] == "random"
        assert "WARNING: enrich_structural_af: 1/1 membership rows without structure " \
               "(method=random, label=negative)" in log


def test_af3_confidence_maps_tokens_and_reorders_chain_pair_iptm():
    with tempfile.TemporaryDirectory() as tmp:
        db, ia, ib, ic = af3_fixture(tmp)
        enrich(tmp, db)
        attrs, (pae, _) = stored(tmp, ia, ib)
        assert attrs["af_source"] == "inferred"
        assert (attrs["iptm"], attrs["ptm"], attrs["ranking_score"]) == (0.81, 0.77, 0.9)
        # Row protein A (P1) is file chain B: the matrix is flipped on both axes.
        np.testing.assert_array_equal(attrs["chain_pair_iptm"], [[0.95, 0.6], [0.7, 0.9]])
        np.testing.assert_array_equal(pae, PAE[np.ix_(ROW_IDX, ROW_IDX)])
        cov = {(r["method"], r["split"], r["label"]): r for r in read_tsv(os.path.join(tmp, "coverage.tsv"))}
        assert cov[("random", "train", "positive")]["with_iptm"] == "1"


def test_missing_pae_is_counted_not_fatal_but_misshapen_pae_is_fatal():
    with tempfile.TemporaryDirectory() as tmp:
        db, ia, ib, _ = afdb_fixture(tmp, write_pae=False)
        log = enrich(tmp, db).stdout
        attrs, pae = stored(tmp, ia, ib)
        assert pae is None and math.isnan(attrs["iptm"])
        assert "1 confidence file(s) missing" in log
    with tempfile.TemporaryDirectory() as tmp:
        db, *_ = afdb_fixture(tmp, pae=PAE[:-1, :-1])
        res = enrich(tmp, db, check=False)
        assert res.returncode != 0 and "does not fit the model" in res.stderr


def test_missing_complex_file_is_reported_as_file_missing():
    with tempfile.TemporaryDirectory() as tmp:
        db, *_ = afdb_fixture(tmp)
        os.remove(os.path.join(tmp, "AF-P2-P1-model_v1.pdb.zst"))
        enrich(tmp, db)
        todo = {(r["uniprot_id_a"], r["uniprot_id_b"]): r for r in read_tsv(os.path.join(tmp, "af3", "all_methods.tsv"))}
        assert todo[("P1", "P2")]["reason"] == "file_missing"
        assert todo[("P1", "P2")]["methods"] == "random,random_ilp"
        assert todo[("P1", "P2")]["splits"] == "test_uniform_1_1,test_uniform_1_10,train"
        assert todo[("P1", "P2")]["n_instance_pairs"] == "1"
