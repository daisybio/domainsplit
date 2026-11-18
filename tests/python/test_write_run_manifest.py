#!/usr/bin/env python3
"""Checks for bin/write_run_manifest.py: hashes, per-split counts, base passthrough."""

import hashlib
import json
import os
import sqlite3
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "bin"))

import write_run_manifest  # noqa: E402


def test_manifest():
    with tempfile.TemporaryDirectory() as d:
        db = os.path.join(d, "m.sqlite3")
        conn = sqlite3.connect(db)
        conn.executescript(
            "CREATE TABLE domain_domain_interaction(id INTEGER PRIMARY KEY, negative INTEGER, source TEXT);"
            "CREATE TABLE ddi_split_membership(ddi_id INTEGER, method TEXT, split TEXT);"
            "INSERT INTO domain_domain_interaction VALUES (1,0,'3did'),(2,1,'sampled_negative'),"
            " (3,1,'sampled_negative'),(4,0,'3did'),"
            " (5,1,'sampled_negative'),(6,1,'Negatome, sampled_negative'),(7,1,'Negatome'),(8,1,'sampledXnegative');"
            # DDI 1 has three instance pairs, DDIs 2 and 3 one each: family ratio 2.0, instance ratio 2/3.
            "INSERT INTO ddi_split_membership VALUES"
            " (1,'random','train'),(1,'random','train'),(1,'random','train'),"
            " (2,'random','train'),(3,'random','train'),"
            " (4,'random_hrni','validation');"
        )
        conn.commit()
        conn.close()

        manual = os.path.join(d, "hippie.tsv")
        with open(manual, "wb") as f:
            f.write(b"a\tb\n")
        base = os.path.join(d, "base.json")
        with open(base, "w") as f:
            json.dump({"params": {"seed": 42}, "urls": {"url_3did": "https://example.invalid/x"}}, f)
        out = os.path.join(d, "run_manifest.json")

        write_run_manifest.main([
            "--base", base, "--db", db, "--out", out, "--no-network",
            "--input", f"hippie_tsv={manual}",
        ])
        m = json.load(open(out))

    assert m["params"]["seed"] == 42
    assert m["manual_inputs"]["hippie_tsv"]["sha256"] == hashlib.sha256(b"a\tb\n").hexdigest()
    assert m["split_counts"]["random"]["train"] == {
        "positives": 1, "negatives": 2, "negatives_per_positive": 2.0,
        "positive_instance_pairs": 3, "negative_instance_pairs": 2, "instance_pairs_negatives_per_positive": 0.6667,
    }
    assert m["sampled_negatives_without_membership"] == 2  # DDIs 5 and 6; 7 is not protected, 8 is no exact token
    assert m["split_counts"]["random_hrni"]["validation"]["negatives_per_positive"] == 0.0
    assert m["url_headers"]["url_3did"] == {"url": "https://example.invalid/x"}


if __name__ == "__main__":
    test_manifest()
    print("ok")
