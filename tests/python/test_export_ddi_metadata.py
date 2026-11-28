#!/usr/bin/env python3
"""Checks for bin/export_ddi_metadata.py against hand-computed values on a
two-residue-per-chain slice.

Slice (row orientation, chain A = domain a):
  A10 CB (0,0,0)    pLDDT 80      B20 CB (5,0,0)   pLDDT 90
  A11 CB (100,0,0)  pLDDT 50      B21 CB (0,7.5,0) pLDDT 70
Contacts (< 8 A): A10-B20 (5), A10-B21 (7.5). PAE rows/cols [A10, A11, B20, B21],
asymmetric: AB A10->B20 = 2, A10->B21 = 12 (on the strict cutoff, excluded),
BA B20->A10 = 20, B21->A10 = 4, everything else 30.
"""

import csv
import gzip
import math
import os
import sqlite3
import sys

import h5py
import numpy as np
import pytest

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import export_ddi_metadata as edm  # noqa: E402
import domainsplit_schema  # noqa: E402

IA, IB = "PF00001_P11111_1_50", "PF00002_P22222_1_50"


def atom(serial, name, resn, chain, resseq, xyz, b):
    x, y, z = xyz
    return (f"ATOM  {serial:5d} {name:<4s} {resn:3s} {chain}{resseq:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{b:6.2f}           {name[0]}")


# (chain, resseq, ca_xyz | None, cb_xyz, plddt)
RESIDUES = [("A", 10, (0, 0, 1), (0, 0, 0), 80), ("A", 11, (100, 0, 1), (100, 0, 0), 50),
            ("B", 20, (5, 0, 1), (5, 0, 0), 90), ("B", 21, (0, 7.5, 1), (0, 7.5, 0), 70)]


def pdb(residues):
    lines, n = [], 0
    for chain, resseq, ca, cb, b in residues:
        if ca is not None:
            n += 1
            lines.append(atom(n, "CA", "ALA", chain, resseq, ca, b))
        n += 1
        lines.append(atom(n, "CB", "ALA", chain, resseq, cb, b))
    return "\n".join(lines + ["END"]) + "\n"


def pae_matrix():
    m = np.full((4, 4), 30.0)
    m[0, 2], m[0, 3], m[2, 0], m[3, 0] = 2, 12, 20, 4
    return m


def d0(n):
    return 1.24 * (max(n, 19) - 15) ** (1 / 3) - 1.8


D0 = d0(2)
TM = lambda p, d: 1 / (1 + (p / d) ** 2)  # noqa: E731
EXPECTED = {
    "min_cb_dist": 5.0,
    "pDockQ": 0.724 / (1 + math.exp(-0.052 * (80 * math.log10(2) - 152.611))) + 0.018,
    "pDockQ2_AB": 1.31034 / (1 + math.exp(-0.075 * (80 * (TM(2, 10) + TM(12, 10)) / 2 - 84.733))) + 0.005,
    "pDockQ2_BA": 1.31034 / (1 + math.exp(-0.075 * ((90 * TM(20, 10) + 70 * TM(4, 10)) / 2 - 84.733))) + 0.005,
    "LIS_AB": 1 - 2 / 12,
    "LIS_BA": 1 - 4 / 12,
    "ipSAE_AB": TM(2, D0),
    "ipSAE_BA": TM(4, D0),
    "ipTM_d0chn_AB": (TM(2, D0) + TM(30, D0)) / 2,
    "ipTM_d0chn_BA": (TM(20, D0) + TM(4, D0)) / 2,
    "numberOfInteractions": 2,
}


def build(tmp_path, residues=RESIDUES, pae=None, res_a=(10, 11), res_b=(20, 21),
          membership=(("m", "train", 1, 1),), with_pae=True):
    db = tmp_path / "db.sqlite3"
    domainsplit_schema.make_db(db)
    conn = sqlite3.connect(db)
    domainsplit_schema.add_instance(conn, "PF00001", "P11111", 1, 50)
    domainsplit_schema.add_instance(conn, "PF00002", "P22222", 1, 50)
    conn.execute("INSERT INTO domain_domain_interaction(id, domain_id_a, domain_id_b, negative) "
                 "SELECT 1, a.id, b.id, 0 FROM domain a, domain b "
                 "WHERE a.pfam_id = 'PF00001' AND b.pfam_id = 'PF00002'")
    for method, split, has, raw in membership:
        conn.execute("INSERT INTO ddi_split_membership(ddi_id, method, split, instance_id_a, instance_id_b, "
                     "has_structure, has_structure_raw) VALUES (1, ?, ?, ?, ?, ?, ?)",
                     (method, split, IA, IB, has, raw))
    conn.commit()
    conn.close()
    h5 = tmp_path / "structures.h5"
    with h5py.File(h5, "w") as f:
        d = f.require_group("PF00001_PF00002").create_dataset(
            f"{IA}__{IB}", data=np.frombuffer(gzip.compress(pdb(residues).encode()), dtype="uint8"))
        d.attrs["af_source"] = "inferred"
        if with_pae:
            p = f.require_group("pae/PF00001_PF00002").create_dataset(
                f"{IA}__{IB}", data=(pae_matrix() if pae is None else pae).astype(np.float32))
            p.attrs["residues_a"] = np.array(res_a, dtype=np.int32)
            p.attrs["residues_b"] = np.array(res_b, dtype=np.int32)
    return db, h5


def run(tmp_path, db, h5):
    out = tmp_path / "out"
    edm.main(["--db", str(db), "--structures_h5", str(h5) if h5 else "", "--outdir", str(out),
              "--coverage", str(tmp_path / "cov.tsv")])
    files = {}
    for method in os.listdir(out):
        with open(out / method / "ddi_metadata.csv") as fh:
            files[method] = list(csv.DictReader(fh))
    with open(tmp_path / "cov.tsv") as fh:
        cov = list(csv.DictReader(fh, delimiter="\t"))
    return files, cov


def check(row, expected):
    for col, want in expected.items():
        got = row[f"af_{col}"]
        if want is None:
            assert got == "", col
        else:
            assert float(got) == pytest.approx(want, abs=1e-3), col


def test_hand_computed_reference(tmp_path):
    files, cov = run(tmp_path, *build(tmp_path))
    (row,) = files["m"]
    assert (row["instance_id_a"], row["instance_id_b"], row["pfam_id_a"], row["uniprot_id_b"], row["source"]) == \
        (IA, IB, "PF00001", "P22222", "inferred")
    check(row, EXPECTED)
    (c,) = cov
    assert (c["source"], c["rows"], c["exported"], c["with_pae"], c["with_contact"]) == ("inferred", "1", "1", "1", "1")


def test_chain_b_written_first_gives_same_scores(tmp_path):
    # slice_pair() keeps the file's chain order, so a swapped slice lists chain B first.
    files, _ = run(tmp_path, *build(tmp_path, residues=RESIDUES[2:] + RESIDUES[:2]))
    check(files["m"][0], EXPECTED)


def test_residue_without_ca_is_dropped_after_alignment(tmp_path):
    # A12 has no CA: it holds a PAE row (zeros, which would change every score
    # if the rows shifted) but no residue in the scoring.
    res = RESIDUES[:2] + [("A", 12, None, (0, 1, 0), 99)] + RESIDUES[2:]
    pae = np.insert(np.insert(pae_matrix(), 2, 0.0, axis=0), 2, 0.0, axis=1)
    files, _ = run(tmp_path, *build(tmp_path, residues=res, pae=pae, res_a=(10, 11, 12)))
    check(files["m"][0], EXPECTED)


def test_contact_cutoff_is_strict(tmp_path):
    # B21 moved to exactly 8 A: A10-B20 is the only contact left.
    res = RESIDUES[:3] + [("B", 21, (0, 8, 1), (0, 8, 0), 70)]
    files, _ = run(tmp_path, *build(tmp_path, residues=res))
    row = files["m"][0]
    assert float(row["af_numberOfInteractions"]) == 1
    assert float(row["af_ipSAE_BA"]) == 0.0  # B20->A10 = 20 is above the PAE cutoff


def test_no_contact_leaves_contact_scores_undefined(tmp_path):
    far = [(c, r, (x + 50 if c == "B" else x, y, z) if ca else None, (cb[0] + 50, cb[1], cb[2]) if c == "B" else cb, b)
           for c, r, ca, cb, b in RESIDUES for x, y, z in [ca]]
    files, cov = run(tmp_path, *build(tmp_path, residues=far))
    check(files["m"][0], {"pDockQ": None, "pDockQ2_AB": None, "pDockQ2_BA": None, "ipSAE_AB": None,
                          "ipSAE_BA": None, "numberOfInteractions": 0, "LIS_AB": 1 - 2 / 12})
    assert cov[0]["with_contact"] == "0"


def test_missing_pae_keeps_only_distance_and_pdockq(tmp_path):
    files, cov = run(tmp_path, *build(tmp_path, with_pae=False))
    check(files["m"][0], {**dict.fromkeys(edm.PAE_COLS), "min_cb_dist": 5.0, "pDockQ": EXPECTED["pDockQ"]})
    assert cov[0]["with_pae"] == "0"


def test_residue_list_mismatch_is_fatal(tmp_path):
    with pytest.raises(ValueError, match="do not match"):
        run(tmp_path, *build(tmp_path, res_a=(10, 12)))


def test_masked_in_one_split_drops_pair_from_that_method_only(tmp_path):
    membership = (("ml", "test_matched_1_1", 1, 1), ("ml", "test_uniform_1_10", 0, 1),
                  ("rnd", "train", 1, 1), ("ext", "train", 0, 0))
    files, cov = run(tmp_path, *build(tmp_path, membership=membership))
    assert files["ml"] == [] and files["ext"] == []
    assert len(files["rnd"]) == 1
    assert {(c["method"], c["split"]): c["exported"] for c in cov} == {
        ("ext", "train"): "0", ("ml", "test_matched_1_1"): "0", ("ml", "test_uniform_1_10"): "0",
        ("rnd", "train"): "1"}


def test_without_structure_writes_header_only(tmp_path):
    db, _ = build(tmp_path, membership=(("m", "train", 0, 0),))
    files, _ = run(tmp_path, db, None)
    assert files == {"m": []}
    with open(tmp_path / "out" / "m" / "ddi_metadata.csv") as fh:
        assert fh.readline().startswith("pfam_id_a,pfam_id_b,uniprot_id_a")


def test_structured_row_without_h5_is_fatal(tmp_path):
    db, _ = build(tmp_path)
    with pytest.raises(SystemExit, match="no --structures_h5"):
        run(tmp_path, db, None)


def test_parser_keeps_hetatm_apart_and_one_altloc():
    ca_a = atom(1, "CA", "ALA", "A", 10, (0, 0, 0), 80)
    cb_b = atom(2, "CB", "ALA", "A", 10, (9, 9, 9), 80)
    cb_a = atom(3, "CB", "ALA", "A", 10, (1, 0, 0), 80)
    ca_a, cb_b, cb_a = ca_a[:16] + "A" + ca_a[17:], cb_b[:16] + "B" + cb_b[17:], cb_a[:16] + "A" + cb_a[17:]
    hoh = "HETATM" + atom(4, "O", "HOH", "A", 10, (5, 5, 5), 0)[6:]
    chains = edm.parse_slice("\n".join([ca_a, cb_b, cb_a, hoh]))
    assert [r[0] for r in chains["A"]] == [10, 10]  # ALA 10 and HOH 10 stay apart
    assert chains["A"][0][2] == (1.0, 0.0, 0.0)     # CB from altloc A, like the CA
