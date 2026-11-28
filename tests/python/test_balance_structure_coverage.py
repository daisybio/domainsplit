#!/usr/bin/env python3
"""Checks for bin/balance_structure_coverage.py on a hand-built master."""

import csv
import json
import os
import sqlite3
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import balance_structure_coverage as bsc  # noqa: E402
import domainsplit_schema  # noqa: E402


def build(tmp_path, splits):
    """splits: {(method, split): [(negative, has_structure), ...]}, one DDI (one row) each."""
    db = tmp_path / "in.sqlite3"
    domainsplit_schema.make_db(db)
    conn = sqlite3.connect(db)
    ddi = 0
    for (method, split), ddis in splits.items():
        for neg, has in ddis:
            ddi += 1
            ia = domainsplit_schema.add_instance(conn, f"PF{ddi:05d}", "P00001", 1, 50)
            ib = domainsplit_schema.add_instance(conn, f"PF{ddi + 500:05d}", "P00002", 1, 50)
            conn.execute("INSERT INTO domain_domain_interaction(id, domain_id_a, domain_id_b, negative) "
                         "SELECT ?, a.id, b.id, ? FROM domain a, domain b WHERE a.pfam_id = ? AND b.pfam_id = ?",
                         (ddi, neg, f"PF{ddi:05d}", f"PF{ddi + 500:05d}"))
            conn.execute("INSERT INTO ddi_split_membership(ddi_id, method, split, instance_id_a, instance_id_b, "
                         "has_structure, z_score) VALUES (?,?,?,?,?,?,?)",
                         (ddi, method, split, ia, ib, has, 1.5 if has else None))
            if has:
                conn.execute("INSERT INTO ddi_interaction_confirmed VALUES (?, ?, 1, 1)", (ddi, method))
    conn.commit()
    conn.close()
    return db


def run(tmp_path, splits, enabled="true", min_coverage="0.3"):
    db = build(tmp_path, splits)
    out = tmp_path / "out.sqlite3"
    report, summary = tmp_path / "r.tsv", tmp_path / "s.json"
    bsc.main(["--db_in", str(db), "--db_out", str(out), "--seed", "42",
              "--enabled", enabled, "--min_coverage", min_coverage,
              "--report", str(report), "--summary", str(summary)])
    with open(report) as fh:
        rep = {(r["split"], r["label"]): r for r in csv.DictReader(fh, delimiter="\t")}
    return sqlite3.connect(out), rep, json.loads(summary.read_text())


BAL = {("m", "train"): [(0, 1)] * 4 + [(1, 1), (1, 1), (1, 0), (1, 0)]}


def test_balanced_masks_higher_label_down(tmp_path):
    conn, rep, summary = run(tmp_path, BAL)
    pos = rep[("train", "positive")]
    assert (pos["covered_before"], pos["covered_after"], pos["masked"], pos["action"]) == ("4", "2", "2", "balanced")
    assert rep[("train", "negative")]["masked"] == "0"
    masked = conn.execute("SELECT ddi_id FROM ddi_split_membership WHERE has_structure_raw = 1 "
                          "AND has_structure = 0").fetchall()
    assert len(masked) == 2
    assert conn.execute("SELECT COUNT(*) FROM ddi_split_membership WHERE has_structure = 0 "
                        "AND z_score IS NOT NULL").fetchone()[0] == 0
    left = {d for (d,) in conn.execute("SELECT ddi_id FROM ddi_interaction_confirmed")}
    assert not left & {d for (d,) in masked} and len(left) == 4
    # Kept = the two lowest hash ranks of the covered positives.
    keys = sorted((bsc.zlib.crc32(f"42:PF{d:05d}:PF{d + 500:05d}".encode()), d) for d in range(1, 5))
    assert {d for (d,) in masked} == {d for _k, d in keys[2:]}
    assert summary["masked_ddis"] == 2 and summary["actions"] == {"m/train": "balanced"}


def test_low_coverage_masks_whole_split(tmp_path):
    splits = {("m", "test_uniform_1_10"): [(0, 1)] * 3 + [(1, 0)] * 3}   # neg coverage 0
    conn, rep, summary = run(tmp_path, splits)
    assert rep[("test_uniform_1_10", "positive")]["action"] == "masked_all"
    assert rep[("test_uniform_1_10", "positive")]["covered_after"] == "0"
    assert conn.execute("SELECT SUM(has_structure), SUM(has_structure_raw) "
                        "FROM ddi_split_membership").fetchone() == (0, 3)


def test_threshold_boundary(tmp_path):
    # 3/10 = 0.3 is not below 0.3: balanced, not masked_all.
    splits = {("m", "train"): [(0, 1)] * 3 + [(0, 0)] * 7 + [(1, 1)] * 3 + [(1, 0)] * 7}
    _conn, rep, _s = run(tmp_path, splits)
    assert rep[("train", "positive")]["action"] == "balanced"
    assert rep[("train", "positive")]["masked"] == "0"


def test_off_reports_only(tmp_path):
    conn, rep, summary = run(tmp_path, BAL, enabled="false")
    assert rep[("train", "positive")]["action"] == "off"
    assert summary["masked_ddis"] == 0
    assert conn.execute("SELECT SUM(has_structure), SUM(has_structure_raw) "
                        "FROM ddi_split_membership").fetchone() == (6, 6)
