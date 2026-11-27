#!/usr/bin/env python3
"""Checks for bin/export_domain_metadata.py on a hand-built master."""

import csv
import os
import sqlite3
import sys

import pytest
from Bio.SeqUtils.ProtParam import ProteinAnalysis

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import export_domain_metadata as edm  # noqa: E402
import domainsplit_schema  # noqa: E402

SEQ = {"P11111": "MKDEARLLGHWYVSTNQCP", "P22222": "MKXAAAAAA", "P33333": "MAZBUOKDE"}


def build(tmp_path, rows, seqs=SEQ):
    db = tmp_path / "db.sqlite3"
    domainsplit_schema.make_db(db)
    conn = sqlite3.connect(db)
    ids = {}
    for name, (pfam, acc, s, e) in {"a": ("PF00001", "P11111", 2, 12), "x": ("PF00002", "P22222", 1, 5),
                                     "z": ("PF00003", "P33333", 1, 9), "c": ("PF00004", "P11111", 13, 19)}.items():
        ids[name] = domainsplit_schema.add_instance(conn, pfam, acc, s, e)
    for acc, seq in seqs.items():
        conn.execute("UPDATE protein SET sequence = ? WHERE uniprot_id = ?", (seq, acc))
    conn.execute("INSERT INTO domain_domain_interaction(id, domain_id_a, domain_id_b, negative) VALUES (1, 1, 2, 0)")
    for method, a, b in rows:
        conn.execute("INSERT INTO ddi_split_membership(ddi_id, method, split, instance_id_a, instance_id_b) "
                     "VALUES (1, ?, 'train', ?, ?)", (method, ids[a], ids[b]))
    conn.commit()
    conn.close()
    return db, ids


def run(tmp_path, db):
    out = tmp_path / "out"
    edm.main(["--db", str(db), "--outdir", str(out), "--report", str(tmp_path / "r.tsv")])
    files = {}
    for m in os.listdir(out):
        with open(out / m / "domain_metadata.csv") as fh:
            files[m] = {r["instance_id"]: r for r in csv.DictReader(fh)}
    with open(tmp_path / "r.tsv") as fh:
        report = {r["metric"]: int(r["value"]) for r in csv.DictReader(fh, delimiter="\t")}
    return files, report


def test_features_per_method_and_counts(tmp_path):
    db, ids = build(tmp_path, [("m1", "a", "x"), ("m1", "z", "a"), ("m2", "c", "a")])
    files, report = run(tmp_path, db)
    assert set(files["m1"]) == {ids["a"], ids["z"]}  # x has an X: skipped
    assert set(files["m2"]) == {ids["a"], ids["c"]}
    row = files["m1"][ids["a"]]
    dom = SEQ["P11111"][1:12]  # 1-based inclusive 2..12
    assert int(row["sequence_length"]) == 11
    assert float(row["isoelectric_point"]) == pytest.approx(ProteinAnalysis(dom).isoelectric_point())
    assert float(row["rel_negative_charge"]) == pytest.approx(2 / 11)  # D, E
    assert float(row["average_flexibility"]) > 0
    assert "sequence" not in row and "isolectric_point" not in row
    # Z/B/U/O substituted (Q/N/C/K), counted once per instance each.
    z = files["m1"][ids["z"]]
    assert float(z["rel_positive_charge"]) == pytest.approx(2 / 9)  # O->K, K
    assert z["average_flexibility"] == ""  # 9 residues: no full flexibility window
    assert report == {"instances": 4, "written": 3, "skipped_x": 1, "substituted_Z": 1,
                      "substituted_B": 1, "substituted_U": 1, "substituted_O": 1}


def test_window_outside_sequence_is_fatal(tmp_path):
    db, _ = build(tmp_path, [("m1", "a", "c")], seqs={**SEQ, "P11111": "MKDEARLLGHWYVST"})
    with pytest.raises(SystemExit, match="outside their protein.sequence"):
        run(tmp_path, db)
