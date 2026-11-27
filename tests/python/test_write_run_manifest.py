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
            "CREATE TABLE domain(id INTEGER PRIMARY KEY, pfam_id TEXT);"
            "INSERT INTO domain VALUES (1,'PF00001'),(2,'PF00002');"
            "CREATE TABLE protein(id INTEGER PRIMARY KEY, uniprot_id TEXT, sequence TEXT);"
            "INSERT INTO protein VALUES (1,'P1','MKV'),(2,'P2',NULL);"
            "CREATE TABLE domain_domain_interaction(id INTEGER PRIMARY KEY, domain_id_a, domain_id_b,"
            " negative INTEGER, source TEXT);"
            "CREATE TABLE ddi_split_membership(ddi_id INTEGER, method TEXT, split TEXT,"
            " instance_id_a TEXT, instance_id_b TEXT);"
            "INSERT INTO domain_domain_interaction VALUES (1,1,2,0,'3did'),(2,1,1,1,'sampled_negative'),"
            " (3,2,2,1,'sampled_negative'),(4,1,2,0,'3did'),"
            " (5,1,2,1,'sampled_negative'),(6,1,2,1,'Negatome, sampled_negative'),(7,1,2,1,'Negatome'),"
            " (8,1,2,1,'sampledXnegative');"
            # DDI 1 has three instance pairs, DDIs 2 and 3 one each: family ratio 2.0, instance ratio 2/3.
            "INSERT INTO ddi_split_membership VALUES"
            " (1,'random','train','PF00001_P1_1_9','PF00002_P2_1_9'),(1,'random','train','PF00001_P3_1_9','PF00002_P2_1_9'),"
            " (1,'random','train','PF00001_P4_1_9','PF00002_P2_1_9'),"
            " (2,'random','train','PF00001_P1_1_9','PF00001_P3_1_9'),(3,'random','train',NULL,NULL),"
            " (4,'random_hrni','validation','PF00001_P1_1_9','PF00002_P2_1_9');"
        )
        conn.commit()
        conn.close()

        manual = os.path.join(d, "hippie.tsv")
        with open(manual, "wb") as f:
            f.write(b"a\tb\n")
        base = os.path.join(d, "base.json")
        with open(base, "w") as f:
            json.dump({"params": {"seed": 42}, "urls": {"url_3did": "https://example.invalid/x",
                                                       "url_pfam2go": "/data/pfam2go"}}, f)
        out = os.path.join(d, "run_manifest.json")
        queue = os.path.join(d, "af3_queue.json")
        with open(queue, "w") as f:
            json.dump({"levels": [{"level": 1, "pairs": 3}]}, f)

        write_run_manifest.main([
            "--base", base, "--db", db, "--out", out, "--no-network",
            "--input", f"hippie_tsv={manual}",
            "--input", f"extra/af3_queue.json={queue}",
            # a local url_* arrives staged under another path; the original value is kept
            "--local", f"url_pfam2go={manual}",
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
    assert m["url_headers"]["url_pfam2go"] == {"url": "/data/pfam2go", "sha256": hashlib.sha256(b"a\tb\n").hexdigest()}
    assert m["membership_fingerprint"]["rows"] == 6
    assert len(m["membership_fingerprint"]["sha256"]) == 64
    assert m["protein_fingerprint"]["rows"] == 2
    # Step summaries are hashed and embedded under their manifest key.
    assert m["af3_queue"] == {"levels": [{"level": 1, "pairs": 3}]}
    assert len(m["manual_inputs"]["af3_queue.json"]["sha256"]) == 64
    assert "structure_coverage_balance" not in m


def test_membership_fingerprint_ignores_surrogate_ids_and_order():
    def db(path, ids, order):
        conn = sqlite3.connect(path)
        conn.executescript(
            "CREATE TABLE domain(id INTEGER PRIMARY KEY, pfam_id TEXT);"
            "CREATE TABLE domain_domain_interaction(id INTEGER PRIMARY KEY, domain_id_a, domain_id_b, negative);"
            "CREATE TABLE ddi_split_membership(ddi_id, method, split, instance_id_a, instance_id_b);")
        conn.executemany("INSERT INTO domain VALUES (?,?)", [(ids[0], "PF1"), (ids[1], "PF2")])
        conn.execute("INSERT INTO domain_domain_interaction VALUES (?,?,?,0)", (ids[2], ids[0], ids[1]))
        rows = [(ids[2], "random", "train", "PF1_A_1_2", "PF2_B_1_2"), (ids[2], "random", "test", "PF1_C_1_2", "PF2_B_1_2")]
        conn.executemany("INSERT INTO ddi_split_membership VALUES (?,?,?,?,?)", rows[::order])
        conn.commit()
        conn.close()
        return write_run_manifest.membership_fingerprint(path)

    with tempfile.TemporaryDirectory() as d:
        a = db(os.path.join(d, "a.sqlite3"), (1, 2, 10), 1)
        b = db(os.path.join(d, "b.sqlite3"), (7, 3, 99), -1)
    assert a == b


if __name__ == "__main__":
    test_manifest()
    print("ok")
