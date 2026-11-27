#!/usr/bin/env python3
"""Checks for bin/fetch_interpro_metadata.py with the InterPro API mocked."""

import csv
import json
import os
import sqlite3
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "bin"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import fetch_interpro_metadata as fim  # noqa: E402
import domainsplit_schema  # noqa: E402

ENTRY = {"metadata": {"literature": {"PUB1": {}, "PUB2": {}}, "counters": {
    "subfamilies": 0, "domain_architectures": 5, "interactions": 0, "matches": 9, "pathways": 0,
    "proteins": 8, "proteomes": 3, "sets": 1, "structural_models": {"alphafold": 6}, "structures": 2, "taxa": 4}}}


def build(tmp_path):
    db = tmp_path / "db.sqlite3"
    domainsplit_schema.make_db(db)
    conn = sqlite3.connect(db)
    a = domainsplit_schema.add_instance(conn, "PF00001", "P11111", 1, 3)
    b = domainsplit_schema.add_instance(conn, "PF00002", "P22222", 1, 3)
    c = domainsplit_schema.add_instance(conn, "PF00003", "P33333", 1, 3)
    conn.execute("INSERT INTO domain_domain_interaction(id, domain_id_a, domain_id_b, negative) VALUES (1, 1, 2, 0)")
    for x, y in ((a, b), (c, a)):
        conn.execute("INSERT INTO ddi_split_membership(ddi_id, method, split, instance_id_a, instance_id_b) "
                     "VALUES (1, 'm', 'train', ?, ?)", (x, y))
    conn.commit()
    conn.close()
    return db


class FakeApi:
    def __init__(self, releases=("111.0",)):
        self.releases = list(releases)
        self.calls = []

    def __call__(self, url, **_kw):
        self.calls.append(url)
        if url == fim.API:
            rel = self.releases.pop(0) if len(self.releases) > 1 else self.releases[0]
            return 200, {"databases": {"interpro": {"version": rel}, "pfam": {"version": "39.0"}}}
        if url.endswith("PF00002/"):
            return 204, None
        if url.endswith("PF00003/"):
            e = json.loads(json.dumps(ENTRY))
            e["metadata"]["literature"] = None
            del e["metadata"]["counters"]["structural_models"]
            return 200, e
        return 200, ENTRY


def run(tmp_path, db, cache, release=""):
    out, summ = tmp_path / "out.tsv", tmp_path / "s.json"
    fim.main(["--db", str(db), "--cache_dir", str(cache), "--release", release,
              "--out", str(out), "--summary", str(summ), "--workers", "2"])
    with open(out) as fh:
        rows = {r["pfam_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    return rows, json.loads(summ.read_text())


def test_fetch_cache_and_missing(tmp_path, monkeypatch):
    db, cache = build(tmp_path), tmp_path / "cache"
    api = FakeApi()
    monkeypatch.setattr(fim, "get", api)
    rows, s = run(tmp_path, db, cache)
    assert rows["PF00001"]["interpro_alphafold"] == "6"
    assert rows["PF00001"]["interpro_literature"] == "2"
    assert rows["PF00002"]["interpro_proteins"] == ""  # 204: missing entry, never 0
    assert rows["PF00003"]["interpro_literature"] == "0"  # null literature: no cited refs
    assert rows["PF00003"]["interpro_alphafold"] == ""  # counter absent: empty
    assert s["families"] == 3 and s["fetched"] == 3 and s["missing_entries"] == 1
    assert s["missing_counter_values"]["alphafold"] == 1
    assert sorted(os.listdir(cache)) == ["interpro-111.0"]

    # Warm run: only the API root is read.
    api.calls.clear()
    rows2, s2 = run(tmp_path, db, cache)
    assert api.calls == [fim.API] and s2["from_cache"] == 3 and rows2 == rows


def test_pinned_release_reads_cache_only(tmp_path, monkeypatch):
    db, cache = build(tmp_path), tmp_path / "cache"
    monkeypatch.setattr(fim, "get", FakeApi())
    rows, _ = run(tmp_path, db, cache)
    api = FakeApi(("112.0",))
    monkeypatch.setattr(fim, "get", api)
    rows2, s = run(tmp_path, db, cache, release="111.0")
    assert rows2 == rows and s["interpro_release"] == "111.0" and s["live_interpro_release"] == "112.0"
    assert api.calls == [fim.API]
    os.remove(cache / "interpro-111.0" / "PF00003.json")
    with pytest.raises(SystemExit, match="missing from"):
        run(tmp_path, db, cache, release="111.0")


def test_pinned_release_without_cache_dir_is_fatal(tmp_path, monkeypatch):
    monkeypatch.setattr(fim, "get", FakeApi(("112.0",)))
    with pytest.raises(SystemExit, match="cache_dir"):
        run(tmp_path, build(tmp_path), "", release="111.0")


def test_release_change_mid_run_is_fatal_and_caches_nothing(tmp_path, monkeypatch):
    monkeypatch.setattr(fim, "get", FakeApi(("111.0", "112.0")))
    with pytest.raises(SystemExit, match="while fetching"):
        run(tmp_path, build(tmp_path), tmp_path / "cache")
    assert not (tmp_path / "cache" / "interpro-111.0").exists()


def test_invalid_cache_record_is_refetched_live_and_fatal_pinned(tmp_path, monkeypatch):
    db, cache = build(tmp_path), tmp_path / "cache"
    monkeypatch.setattr(fim, "get", FakeApi())
    run(tmp_path, db, cache)
    bad = cache / "interpro-111.0" / "PF00001.json"
    bad.write_text('{"status": "ok", "counters": {"proteins": "8"}')  # truncated
    api = FakeApi()
    monkeypatch.setattr(fim, "get", api)
    rows, s = run(tmp_path, db, cache)
    assert s["fetched"] == 1 and rows["PF00001"]["interpro_proteins"] == "8"
    assert json.loads(bad.read_text())["counters"]["proteins"] == 8
    bad.write_text('{"status": "ok", "counters": {"proteins": "8"}}')  # wrong type
    monkeypatch.setattr(fim, "get", FakeApi(("112.0",)))
    with pytest.raises(SystemExit, match="invalid in"):
        run(tmp_path, db, cache, release="111.0")


def test_bad_release_pin_is_fatal(tmp_path, monkeypatch):
    monkeypatch.setattr(fim, "get", FakeApi())
    with pytest.raises(SystemExit, match="not a release number"):
        run(tmp_path, build(tmp_path), tmp_path / "cache", release="../x")


class FakeResponse:
    def __init__(self, status, body):
        self.status, self.body = status, body

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def read(self):
        return self.body


def test_get_204_is_missing_but_empty_200_is_retried(monkeypatch):
    responses = [FakeResponse(204, b""), FakeResponse(200, b""), FakeResponse(200, b'{"metadata": 1}')]
    monkeypatch.setattr(fim.urllib.request, "urlopen", lambda *a, **k: responses.pop(0))
    monkeypatch.setattr(fim.time, "sleep", lambda s: None)
    assert fim.get("https://x/") == (204, None)
    assert fim.get("https://x/") == (200, {"metadata": 1})
