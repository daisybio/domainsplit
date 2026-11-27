#!/usr/bin/env python3
"""Checks for bin/report_split_qc.py on a hand-built master."""

import csv
import os
import sqlite3
import sys

import h5py
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "bin"))

import report_split_qc  # noqa: E402


def inst(pfam, acc, s=1, e=50):
    return f"{pfam}_{acc}_{s}_{e}"


def make_db(path, membership):
    """membership: (ddi_id, method, split, a, b); DDIs 1-3 positive, 4-6 negative."""
    conn = sqlite3.connect(path)
    conn.executescript(
        "CREATE TABLE domain_domain_interaction(id INTEGER PRIMARY KEY, negative INTEGER, source TEXT);"
        "CREATE TABLE ddi_split_membership(ddi_id, method, split, instance_id_a, instance_id_b,"
        " has_structure INTEGER NOT NULL DEFAULT 0, z_score REAL);"
        "INSERT INTO domain_domain_interaction VALUES (1,0,'3did'),(2,0,'3did'),(3,0,'PPIDM_Gold'),"
        " (4,1,'sampled_negative'),(5,1,'sampled_negative'),(6,1,'INVERSE_PPIDM_Gold');"
    )
    conn.executemany("INSERT INTO ddi_split_membership(ddi_id, method, split, instance_id_a, instance_id_b)"
                     " VALUES (?,?,?,?,?)", membership)
    conn.commit()
    conn.close()


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


A, B, C, D = inst("PF1", "P1"), inst("PF2", "P2"), inst("PF3", "P3"), inst("PF4", "P4")
E = inst("PF5", "P5")

CLEAN = [
    (1, "ml", "train", A, B), (4, "ml", "train", A, A),
    (2, "ml", "test_matched_1_1", C, D), (5, "ml", "test_matched_1_1", C, C),
    (2, "ml", "test_uniform_1_10", C, D), (5, "ml", "test_uniform_1_10", D, D),
    (1, "ext", "train", A, B), (2, "ext", "train", A, C),
    (3, "ext", "test", A, D), (6, "ext", "test", B, E),
]


def run(tmp_path, membership, *extra):
    db = str(tmp_path / "m.sqlite3")
    make_db(db, membership)
    return report_split_qc.main(["--db", db, "--ilp-methods", "ml", "--external-methods", "ext",
                                 "--outdir", str(tmp_path), *extra])


def test_clean_run(tmp_path):
    assert run(tmp_path, CLEAN) == 0
    assert read(tmp_path / "family_leaks.tsv") == []

    ov = {(r["method"], r["split_a"], r["split_b"], r["unit"]): r for r in read(tmp_path / "split_overlap.tsv")}
    r = ov["ml", "test_matched_1_1", "test_uniform_1_10", "family"]
    assert (r["n_a"], r["n_b"], r["n_shared"]) == ("2", "2", "2")
    assert ov["ml", "test_matched_1_1", "train", "instance"]["n_shared"] == "0"

    var = {(r["label"], r["granularity"]): r for r in read(tmp_path / "test_variant_overlap.tsv")}
    assert var["positive", "ddi_row"]["n_shared"] == "1"
    assert var["negative", "ddi_row"]["n_shared"] == "0"
    assert var["negative", "family_pair"]["n_1_1"] == "1"
    assert {r["method"] for r in read(tmp_path / "test_variant_overlap.tsv")} == {"ml"}

    prof = {(r["label"], r["source"]): r for r in read(tmp_path / "external_test_profile.tsv")}
    pos = prof["positive", "PPIDM_Gold"]
    # PF1 has train-positive degree 2 (PF2, PF3), PF4 none.
    assert (pos["families"], pos["families_in_train"]) == ("2", "1")
    assert (pos["train_pos_degree_min"], pos["train_pos_degree_max"]) == ("0", "2")
    assert pos["pair_min_train_pos_degree_median"] == "0"
    neg = prof["negative", "__total__"]
    assert neg["families_shared_with_other_label"] == "0"
    assert neg["train_pos_degree_max"] == "1"


def test_leak_fails(tmp_path):
    leaky = CLEAN + [(4, "ml", "validation", A, C)]
    assert run(tmp_path, leaky) == 1
    leaks = read(tmp_path / "family_leaks.tsv")
    assert [(r["pfam_id"], r["splits"]) for r in leaks] == [
        ("PF1", "train,validation"),
        ("PF3", "test_matched_1_1,test_uniform_1_10,validation"),
    ]


def test_external_test_split_not_checked(tmp_path):
    # ext's test shares PF1 with its train: reported, but no guarantee covers it.
    db = str(tmp_path / "m.sqlite3")
    make_db(db, CLEAN)
    assert report_split_qc.main(["--db", db, "--ilp-methods", "ml,ext", "--external-methods", "ext",
                                 "--outdir", str(tmp_path)]) == 0


def test_missing_method_fails(tmp_path):
    db = str(tmp_path / "m.sqlite3")
    make_db(db, CLEAN)
    try:
        report_split_qc.main(["--db", db, "--ilp-methods", "nope", "--outdir", str(tmp_path)])
    except SystemExit as e:
        assert "nope" in str(e)
    else:
        raise AssertionError("an ILP method with no rows must fail")


def test_esm3_coverage(tmp_path):
    h5 = tmp_path / "esm3_structure_domain_embeddings.h5"
    with h5py.File(h5, "w") as f:
        for i in (A, B, C):
            f.create_dataset(f"{i.split('_')[0]}/{i}", data=np.zeros(4))
    report = tmp_path / "afdb_structures_report.tsv"
    report.write_text("uniprot_id\tstatus\tlength\tplaced\tmismatch_rate\tcoverage\tdetail\n"
                      "P1\tok\t\t\t\t\t\nP4\tok\t\t\t\t\t\nP5\tnot_in_afdb\t\t\t\t\t\n")
    assert run(tmp_path, CLEAN, "--esm3-h5", str(h5), "--afdb-report", str(report)) == 0

    cov = {(r["method"], r["split"], r["label"]): r for r in read(tmp_path / "esm3_structure_coverage.tsv")}
    assert (cov["ml", "train", "positive"]["both_in_h5"], cov["ml", "train", "positive"]["instance_pairs"]) == ("1", "1")
    assert cov["ml", "test_uniform_1_10", "negative"]["both_in_h5"] == "0"

    miss = {r["instance_id"]: r for r in read(tmp_path / "esm3_structure_missing_instances.tsv")}
    assert set(miss) == {D, E}
    assert miss[D]["reason"] == "below_min_coverage"
    assert miss[D]["n_membership_rows"] == "4"
    assert miss[E]["reason"] == "not_in_afdb"
    acc = read(tmp_path / "esm3_structure_missing_accessions.tsv")
    assert [(r["accession"], r["reason"], r["n_instances"]) for r in acc] == [
        ("P4", "below_min_coverage", "1"), ("P5", "not_in_afdb", "1")]


def test_z_score_summary(tmp_path):
    db = str(tmp_path / "m.sqlite3")
    make_db(db, CLEAN)
    conn = sqlite3.connect(db)
    # ml train positive scored (2.0); ml train negative structured but unscored.
    conn.execute("UPDATE ddi_split_membership SET has_structure = 1, z_score = 2.0 WHERE method = 'ml' AND ddi_id = 1")
    conn.execute("UPDATE ddi_split_membership SET has_structure = 1 WHERE method = 'ml' AND ddi_id = 4")
    conn.commit()
    conn.close()
    assert report_split_qc.main(["--db", db, "--ilp-methods", "ml", "--external-methods", "ext",
                                 "--outdir", str(tmp_path)]) == 0
    z = {(r["method"], r["split"], r["label"]): r for r in read(tmp_path / "z_score_summary.tsv")}
    pos, neg = z["ml", "train", "positive"], z["ml", "train", "negative"]
    assert (pos["rows"], pos["with_structure"], pos["scored"], pos["z_median"]) == ("1", "1", "1", "2.0")
    assert (neg["scored"], neg["structure_no_score"], neg["z_median"]) == ("0", "1", "")
    assert pos["with_structure_raw"] == ""  # this master has no has_structure_raw column
    assert z["ext", "test", "negative"]["with_structure"] == "0"
