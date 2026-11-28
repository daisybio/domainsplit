#!/usr/bin/env python3
"""Checks for bin/frozen_splits.py: export layout/manifest, verify gate."""

import json
import os
import sys

import sqlite3

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "bin"))

import frozen_splits  # noqa: E402


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return str(path)


def master(path, extra_row=False, sequence="MKV"):
    conn = sqlite3.connect(path)
    conn.executescript(
        "CREATE TABLE domain(id INTEGER PRIMARY KEY, pfam_id);"
        "CREATE TABLE domain_domain_interaction(id INTEGER PRIMARY KEY, domain_id_a, domain_id_b, negative);"
        "CREATE TABLE ddi_split_membership(ddi_id, method, split, instance_id_a, instance_id_b);"
        "CREATE TABLE protein(id INTEGER PRIMARY KEY, uniprot_id, sequence);"
        "INSERT INTO domain VALUES (1,'PF1'),(2,'PF2');"
        "INSERT INTO domain_domain_interaction VALUES (1,1,2,0);"
        "INSERT INTO ddi_split_membership VALUES (1,'random','train','PF1_P1_1_2','PF2_P2_1_2');")
    if extra_row:
        conn.execute("INSERT INTO ddi_split_membership VALUES (1,'external_test','test','PF1_P3_1_2','PF2_P2_1_2')")
    conn.execute("INSERT INTO protein VALUES (1,'P1',?)", (sequence,))
    conn.commit()
    conn.close()
    return str(path)


@pytest.fixture
def bundle(tmp_path):
    src = tmp_path / "src"
    ddis = write(src / "split_ddis.csv", "a,b\nPF1,PF2\nPF3,PF4\n")
    fetch = [write(src / n, n) for n in ("sequences.fasta", "instances.tsv")]
    fam = [f"random:train:{write(src / 'f1' / 'train.csv', 'x')}",
           f"random:validation:{write(src / 'f2' / 'val.csv', 'y')}"]
    inst = [f"random:train:{write(src / 'i1' / 'train.csv', 'xi')}",
            f"random:validation:{write(src / 'i2' / 'val.csv', 'yi')}"]
    out = tmp_path / "frozen_splits"
    argv = ["export", "--split_ddis", ddis, "--params", json.dumps({"seed": "42"}), "--out", str(out),
            "--db", master(tmp_path / "m.sqlite3")]
    for f in fetch:
        argv += ["--fetch", f]
    for f in fam:
        argv += ["--family", f]
    for f in inst:
        argv += ["--instance", f]
    frozen_splits.main(argv)
    return out, ddis


def test_export_layout(bundle):
    out, _ = bundle
    m = json.load(open(out / "bundle_manifest.json"))
    assert m["family"] == [{"method": "random", "split": "train", "path": "splits/family/random/train.csv"},
                           {"method": "random", "split": "validation", "path": "splits/family/random/validation.csv"}]
    assert open(out / "splits/instance/random/validation.csv").read() == "yi"
    assert m["fetch"]["instances.tsv"] == "fetch/instances.tsv"
    assert m["params"] == {"seed": "42"}
    assert set(m["files"]) == {"fetch/sequences.fasta", "fetch/instances.tsv", "split_ddis/split_ddis.csv",
                               "splits/family/random/train.csv", "splits/family/random/validation.csv",
                               "splits/instance/random/train.csv", "splits/instance/random/validation.csv"}


def test_verify_accepts_reordered_population(bundle, tmp_path, capsys):
    out, _ = bundle
    now = write(tmp_path / "now.csv", "a,b\nPF3,PF4\nPF1,PF2\n")
    frozen_splits.main(["verify", "--bundle", str(out), "--split_ddis", now])
    assert "bundle verified" in capsys.readouterr().out


def test_verify_rejects_changed_population(bundle, tmp_path):
    out, _ = bundle
    now = write(tmp_path / "now.csv", "a,b\nPF1,PF2\nPF5,PF6\n")
    with pytest.raises(SystemExit, match="1 row\\(s\\) only frozen .* 1 only in this run"):
        frozen_splits.main(["verify", "--bundle", str(out), "--split_ddis", now])


def test_verify_integrity_only(bundle, capsys):
    out, _ = bundle
    frozen_splits.main(["verify", "--bundle", str(out)])
    assert "integrity verified" in capsys.readouterr().out


def test_verify_rejects_tampered_file(bundle):
    out, ddis = bundle
    write(out / "splits/family/random/train.csv", "changed")
    with pytest.raises(SystemExit, match="missing or changed"):
        frozen_splits.main(["verify", "--bundle", str(out), "--split_ddis", ddis])


def test_export_rejects_mismatched_levels(tmp_path):
    ddis = write(tmp_path / "d.csv", "a\n")
    f = write(tmp_path / "f.csv", "x")
    with pytest.raises(SystemExit, match="different"):
        frozen_splits.main(["export", "--split_ddis", ddis, "--params", "{}", "--out", str(tmp_path / "o"),
                            "--db", master(tmp_path / "m.sqlite3"), "--family", f"random:train:{f}"])


def test_check_db(bundle, tmp_path, capsys):
    out, _ = bundle
    frozen_splits.main(["check_db", "--bundle", str(out), "--db", master(tmp_path / "same.sqlite3")])
    assert "matches the frozen run" in capsys.readouterr().out
    with pytest.raises(SystemExit, match="membership rows 1 -> 2"):
        frozen_splits.main(["check_db", "--bundle", str(out), "--db", master(tmp_path / "x.sqlite3", extra_row=True)])
    with pytest.raises(SystemExit, match="protein rows"):
        frozen_splits.main(["check_db", "--bundle", str(out), "--db", master(tmp_path / "y.sqlite3", sequence="MKA")])


def test_verify_rejects_unhashed_reference(bundle):
    out, ddis = bundle
    m = json.load(open(out / "bundle_manifest.json"))
    del m["files"]["splits/instance/random/train.csv"]
    json.dump(m, open(out / "bundle_manifest.json", "w"))
    with pytest.raises(SystemExit, match="without a hash"):
        frozen_splits.main(["verify", "--bundle", str(out), "--split_ddis", ddis])
