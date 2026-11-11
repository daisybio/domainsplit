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
            "CREATE TABLE domain_domain_interaction(id INTEGER PRIMARY KEY, negative INTEGER);"
            "CREATE TABLE ddi_split_membership(ddi_id INTEGER, method TEXT, split TEXT);"
            "INSERT INTO domain_domain_interaction VALUES (1,0),(2,1),(3,1),(4,0);"
            "INSERT INTO ddi_split_membership VALUES"
            " (1,'random','train'),(2,'random','train'),(3,'random','train'),"
            " (4,'random_hcni','validation');"
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
    assert m["split_counts"]["random"]["train"] == {"positives": 1, "negatives": 2, "negatives_per_positive": 2.0}
    assert m["split_counts"]["random_hcni"]["validation"]["negatives_per_positive"] == 0.0
    assert m["url_headers"]["url_3did"] == {"url": "https://example.invalid/x"}


if __name__ == "__main__":
    test_manifest()
    print("ok")
