#!/usr/bin/env python3
"""Checks for bin/build_af_metadata.py: sources, precedence, dedup, fingerprint."""

import csv
import os
import sys
import time

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "bin"))

import build_af_metadata  # noqa: E402
from enrich_structural_af import build_meta_index  # noqa: E402


def touch(path, text="x"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return path


def af3_job(root, name):
    d = os.path.join(root, "batch0001", name)
    touch(os.path.join(d, f"{name}_model.cif"))
    touch(os.path.join(d, f"{name}_confidences.json"))
    os.makedirs(os.path.join(d, "seed-1_sample-0"), exist_ok=True)
    return os.path.join(d, f"{name}_model.cif")


def afdb_model(root, pair, mid):
    d = os.path.join(root, pair)
    touch(os.path.join(d, f"{mid}-predicted_aligned_error_v1.json.zst"))
    return touch(os.path.join(d, f"{mid}-model_v1.pdb.zst")), d


def run(tmp_path, *argv):
    out = {k: str(tmp_path / f"{k}.out") for k in ("out", "inventory", "report", "digest")}
    build_af_metadata.main([*argv, *[a for k, v in out.items() for a in (f"--{k}", v)]])
    rows = list(csv.DictReader(open(out["out"])))
    inv = list(csv.DictReader(open(out["inventory"]), delimiter="\t"))
    rep = {r["source"]: r for r in csv.DictReader(open(out["report"]), delimiter="\t")}
    return out, rows, inv, rep


@pytest.fixture
def layout(tmp_path):
    dl_root, inf_root = tmp_path / "afdb", tmp_path / "joint"
    m1, d1 = afdb_model(dl_root, "P11111_P22222", "AF-1")
    _m2, d2 = afdb_model(dl_root, "P33333_P44444", "AF-2")
    os.remove(_m2)                                        # downloaded model missing
    manifest = touch(str(tmp_path / "manifest.csv"),
                     "uniprot_id_a,uniprot_id_b,model_entity_id,local_tar_name,path\n"
                     f"P11111,P22222,AF-1,t.tar,{d1}\n"
                     f"P11111,P22222,AF-1,t.tar,{d1}\n"                # exact repeat
                     f"P33333,P44444,AF-2,t.tar,{d2}\n")
    af3_job(inf_root, "Pos_P22222_P11111")                        # reverse order of a downloaded pair
    inf34 = af3_job(inf_root, "Pos_P33333_P44444")                # usable although downloaded is missing
    af3_job(inf_root, "Pos_P55555_P66666-2")                      # isoform accession
    os.makedirs(inf_root / "batch0001" / "Pos_P77777_P88888")      # job without a model
    return dict(manifest=manifest, inf_root=str(inf_root), m1=m1, inf34=inf34)


def test_precedence_dedup_and_missing(tmp_path, layout):
    _out, rows, inv, rep = run(tmp_path, "--manifest", layout["manifest"], "--inferred_dir", layout["inf_root"])
    got = {(r["uniprot_id_a"], r["uniprot_id_b"]): (r["source"], r["path"]) for r in rows}
    assert got == {
        ("P11111", "P22222"): ("downloaded", layout["m1"]),       # downloaded wins over Pos_P22222_P11111
        ("P33333", "P44444"): ("inferred", layout["inf34"]),      # missing downloaded does not shadow
        ("P55555", "P66666-2"): ("inferred", got[("P55555", "P66666-2")][1]),
    }
    assert rep["downloaded"] == {"source": "downloaded", "candidates": "3", "model_missing": "1",
                                 "exact_duplicate": "1", "shadowed": "0", "selected": "1"}
    assert rep["inferred"]["shadowed"] == "1" and rep["inferred"]["model_missing"] == "1"
    assert {r["status"] for r in inv} == {"selected", "shadowed", "model_missing", "exact_duplicate"}
    # The job's seed subdirectory was not mistaken for a job.
    assert not any("seed-1" in r["path"] for r in inv)


def test_fingerprint_and_sidecars(tmp_path, layout):
    _out, rows, _inv, _rep = run(tmp_path, "--manifest", layout["manifest"], "--inferred_dir", layout["inf_root"])
    by = {r["uniprot_id_a"]: r for r in rows}
    assert int(by["P11111"]["pae_size"]) == 1 and int(by["P11111"]["summary_size"]) == -1   # AFDB: no summary file
    assert int(by["P33333"]["pae_size"]) == 1 and int(by["P33333"]["summary_size"]) == -1   # AF3 summary absent here


def test_sidecar_change_changes_digest(tmp_path, layout):
    out, *_ = run(tmp_path, "--inferred_dir", layout["inf_root"])
    d1 = open(out["digest"]).read()
    time.sleep(0.01)
    touch(layout["inf34"].replace("_model.cif", "_confidences.json"), "changed")
    out, *_ = run(tmp_path, "--inferred_dir", layout["inf_root"])
    assert open(out["digest"]).read() != d1


def test_same_source_conflict_is_fatal(tmp_path, layout):
    other = tmp_path / "joint2"
    af3_job(str(other), "Pos_P33333_P44444")
    with pytest.raises(SystemExit, match="two different inferred models"):
        run(tmp_path, "--inferred_dir", layout["inf_root"], "--inferred_dir", str(other))


def test_relative_af_metadata_and_reader_compat(tmp_path, layout):
    rel = os.path.relpath(layout["inf34"], tmp_path)
    meta = touch(str(tmp_path / "af.csv"), f"uniprot_id_a,uniprot_id_b,path,source\nP33333,P44444,{rel},inferred\n")
    out, rows, _inv, _rep = run(tmp_path, "--af_metadata", meta)
    assert rows[0]["path"] == layout["inf34"]               # absolute after resolution
    index = build_meta_index(out["out"])                    # extra columns are ignored
    assert index[frozenset({"P33333", "P44444"})] == (layout["inf34"], "P33333", "inferred")


def test_job_name_needs_accessions_and_root_can_be_a_job(tmp_path):
    root = tmp_path / "batch_2026_11"              # three tokens, but not accessions: not a job
    model = af3_job(str(root), "Pos_P12345_Q67890")
    _o, rows, _i, _r = run(tmp_path, "--inferred_dir", str(root))
    assert [(r["uniprot_id_a"], r["path"]) for r in rows] == [("P12345", model)]
    _o, rows, _i, _r = run(tmp_path, "--inferred_dir", os.path.dirname(model))
    assert [r["uniprot_id_b"] for r in rows] == ["Q67890"]


def test_nothing_usable_fails(tmp_path):
    empty = tmp_path / "e"
    empty.mkdir()
    with pytest.raises(SystemExit, match="no usable AF model"):
        run(tmp_path, "--inferred_dir", str(empty))
