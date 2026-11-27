#!/usr/bin/env python3
"""Checks for bin/build_af3_queue.py on a hand-built master."""

import csv
import json
import os
import sqlite3
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import build_af3_queue  # noqa: E402
import domainsplit_schema  # noqa: E402

# (ddi_id, negative, method, split, uniprot_a, uniprot_b, has_structure)
ROWS = [
    (1, 0, "minimal_leakage", "train", "P00002", "P00001", 0),       # level 1, stored B/A
    (2, 1, "random", "train", "P00001", "P00002", 0),                # same pair, later level
    (3, 1, "minimal_leakage", "test_uniform_1_10", "P00003", "P00004", 0),  # level 8 only
    (4, 0, "random", "validation", "P00003", "P00004", 0),           # level 7 first
    (5, 0, "minimal_leakage", "validation", "P00005", "P00006", 1),  # sliced model
    (6, 1, "minimal_leakage", "test_matched_1_1", "P00007", "P00008", 0),  # unusable model
    (7, 0, "random_ilp", "train", "P00009", "P00010", 0),            # method not selected
]


def build(tmp_path, rows=ROWS):
    db = tmp_path / "master.sqlite3"
    domainsplit_schema.make_db(db)
    conn = sqlite3.connect(db)
    for ddi, neg, method, split, ua, ub, has in rows:
        ia = domainsplit_schema.add_instance(conn, f"PF{ddi:05d}", ua, 1, 50)
        ib = domainsplit_schema.add_instance(conn, f"PF{ddi + 100:05d}", ub, 1, 50)
        conn.execute("INSERT OR IGNORE INTO domain_domain_interaction(id, domain_id_a, domain_id_b, negative) "
                     "SELECT ?, a.id, b.id, ? FROM domain a, domain b WHERE a.pfam_id = ? AND b.pfam_id = ?",
                     (ddi, neg, f"PF{ddi:05d}", f"PF{ddi + 100:05d}"))
        conn.execute("INSERT INTO ddi_split_membership(ddi_id, method, split, instance_id_a, instance_id_b, "
                     "has_structure) VALUES (?,?,?,?,?,?)", (ddi, method, split, ia, ib, has))
    conn.commit()
    conn.close()
    meta = tmp_path / "af_metadata.csv"
    meta.write_text("uniprot_id_a,uniprot_id_b,path,source\n"
                    "P00006,P00005,/m/a.cif,inferred\nP00007,P00008,/m/b.pdb.zst,downloaded\n")
    gap = tmp_path / "all_methods.tsv"
    gap.write_text("uniprot_id_a\tuniprot_id_b\treason\tn_instance_pairs\tlabels\tsplits\tmethods\n"
                   "P00007\tP00008\tslice_failed\t1\tnegative\ttest_matched_1_1\tminimal_leakage\n")
    return db, meta, gap


def run(tmp_path, methods="minimal_leakage,random", rows=ROWS):
    db, meta, gap = build(tmp_path, rows)
    out = tmp_path / "out"
    build_af3_queue.main(["--db", str(db), "--af_metadata", str(meta), "--gap_report", str(gap),
                          "--methods", methods, "--seed", "42", "--out_dir", str(out),
                          "--counts", str(out / "af3_queue.json")])
    return out


def data_rows(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_queue_levels_and_dedup(tmp_path):
    out = run(tmp_path)
    rows = data_rows(out / "queue.tsv")
    got = {(r["uniprot_a"], r["uniprot_b"]): r for r in rows}
    # Each missing pair once, at its first level; modelled pairs never.
    assert set(got) == {("P00001", "P00002"), ("P00003", "P00004")}
    assert got[("P00001", "P00002")]["level"] == "1"
    assert got[("P00001", "P00002")]["labels"] == "pos"          # random/train neg is a later level
    assert got[("P00001", "P00002")]["splits"] == "minimal_leakage/train"
    assert got[("P00003", "P00004")]["level"] == "7"
    assert got[("P00003", "P00004")]["level_name"] == "random"


def test_queue_plain_tsv_and_summary(tmp_path):
    out = run(tmp_path)
    lines = (out / "queue.tsv").read_text().splitlines()
    assert lines[0].split("\t") == build_af3_queue.COLUMNS
    assert not any(l.startswith("#") for l in lines)               # data only
    assert [l.split("\t")[0] for l in lines[1:]] == ["1", "7"]     # level order
    summary = (out / "queue_summary.txt").read_text().splitlines()
    assert summary[0] == "minimal_leakage (1 pos, 0 neg)"
    assert summary[1].startswith("+ ") and len(summary) == 3      # minimal_leakage, random, test_uniform_1_10
    counts = json.loads((out / "af3_queue.json").read_text())
    assert [c["pairs"] for c in counts["levels"]] == [1, 0, 0, 0, 0, 0, 1, 0]
    assert counts["unusable_models"] == 1


def test_unusable_models(tmp_path):
    out = run(tmp_path)
    rows = data_rows(out / "unusable_models.tsv")
    assert len(rows) == 1
    r = rows[0]
    assert (r["uniprot_a"], r["uniprot_b"], r["reason"], r["source"], r["path"]) == \
        ("P00007", "P00008", "slice_failed", "downloaded", "/m/b.pdb.zst")
    assert r["labels"] == "neg" and r["level"] == "1"


def test_unknown_split_fails(tmp_path):
    rows = ROWS + [(8, 0, "minimal_leakage", "holdout", "P00011", "P00012", 0)]
    with pytest.raises(SystemExit, match="minimal_leakage/holdout"):
        run(tmp_path, rows=rows)


def test_order_is_seeded_hash():
    pairs = [(f"P{i:05d}", f"Q{i:05d}") for i in range(50)]
    a = sorted(pairs, key=lambda p: build_af3_queue.rank("42", *p))
    b = sorted(pairs, key=lambda p: build_af3_queue.rank("43", *p))
    assert a != sorted(pairs) and a != b
