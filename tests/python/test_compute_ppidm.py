#!/usr/bin/env python3
"""Checks for bin/compute_ppidm.py on a toy network of five sources.

Planted:
  * PF00010 x PF00020 co-occur across proteins in s1-s3 (hi_conf) -> PPIDM Gold;
  * PF00600 x PF00601 only in s1 -> PPIDM Silver;
  * PF00100 and PF00200 sit on many PPIs but never on the two sides of one, in
    s1-s4 -> INVERSE_PPIDM Gold (s5 is empty and still counts in S);
  * PF00900 x PF00901: two rare domains that never co-occur, but no source has
    the power to call that depletion -> in neither output;
  * PF00070 + PF00080 on one protein: same-protein co-occurrence, never K;
  * PF00500 homo PPI and homo 3did row: never emitted.

Run directly (`python3 tests/python/test_compute_ppidm.py`) or via pytest.
"""

import csv
import gzip
import json
import math
import os
import sqlite3
import subprocess
import sys
import tempfile

import numpy as np
from scipy.stats import hypergeom

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

SCRIPT = os.path.join(BIN, "compute_ppidm.py")

import compute_ppidm as cp  # noqa: E402
from ddi_db_utils import ensure_domains, insert_ddis, merge_external_ddis  # noqa: E402
from domainsplit_schema import make_db  # noqa: E402

ENV = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))
SOURCES = ["s1", "s2", "s3", "s4", "s5"]
ALPHA = "1.0"  # alpha_low = 1 / C(D, 2): depletion is reachable at toy scale


def toy_proteins():
    prot = {}
    for i in range(12):
        prot[f"X{i}"] = ["PF00100"]
        prot[f"Y{i}"] = ["PF00200"]
        prot[f"F{i}"] = [f"PF0030{i % 4 + 1}"]
    for i in range(4):
        prot[f"E{i}"], prot[f"G{i}"] = ["PF00010"], ["PF00020"]
        prot[f"P{i}"], prot[f"Q{i}"] = ["PF00600"], ["PF00601"]
    for i in range(3):
        prot[f"A{i}"], prot[f"B{i}"] = ["PF00030"], ["PF00040"]
        prot[f"C{i}"], prot[f"D{i}"] = ["PF00050"], ["PF00060"]
    prot.update(R1=["PF00900"], R2=["PF00901"], M1=["PF00070", "PF00080"],
                H1=["PF00500"], H2=["PF00500"])
    return prot


def toy_edges():
    """{source: [(acc_a, acc_b, hi_conf)]}"""
    core = []
    for i in range(12):
        core += [(f"X{i}", f"F{i % 12}"), (f"X{i}", f"F{(i + 5) % 12}"),
                 (f"Y{i}", f"F{(i + 2) % 12}"), (f"Y{i}", f"F{(i + 7) % 12}")]
    core += [(f"A{i}", f"B{i}") for i in range(3)] + [(f"C{i}", f"D{i}") for i in range(3)]
    core += [("R1", "F0"), ("R2", "F1"), ("M1", "F2"), ("H1", "H2"), ("X0", "X1")]
    enriched = [(f"E{i}", f"G{i}") for i in range(4)]
    edges = {s: [(a, b, 1) for a, b in core + enriched] for s in ("s1", "s2", "s3")}
    edges["s1"] += [(f"P{i}", f"Q{i}", 1) for i in range(4)]
    edges["s4"] = [(a, b, 1) for a, b in core] + [(a, b, 0) for a, b in enriched]
    edges["s5"] = []
    return edges


def write_inputs(tmp):
    with gzip.open(os.path.join(tmp, "protein_pfam.tsv.gz"), "wt") as fh:
        fh.write("acc\tpfams\n")
        for acc, pf in sorted(toy_proteins().items()):
            fh.write(f"{acc}\t{';'.join(sorted(pf))}\n")
    for s, rows in toy_edges().items():
        with gzip.open(os.path.join(tmp, f"source_{s}.tsv.gz"), "wt") as fh:
            fh.write("acc_a\tacc_b\thi_conf\n")
            for a, b, hi in rows:
                a, b = sorted((a, b))
                fh.write(f"{a}\t{b}\t{hi}\n")
    db = os.path.join(tmp, "gold.sqlite3")
    make_db(db)
    conn = sqlite3.connect(db)
    gold = [("PF00030", "PF00040"), ("PF00050", "PF00060"), ("PF00500", "PF00500"),
            ("PF00030", "PF09999")]  # last one: domain outside the universe
    ensure_domains(conn, {p for pair in gold for p in pair})
    insert_ddis(conn, gold, negative=False, source="3did")
    merge_external_ddis(conn, [("PF00100", "PF00200", 1, "negatome")])  # not gold: negative
    conn.commit()
    conn.close()


def run(tmp, out, *extra, sources=SOURCES):
    os.makedirs(out, exist_ok=True)
    args = [sys.executable, SCRIPT, "--protein-pfam", os.path.join(tmp, "protein_pfam.tsv.gz")]
    for s in sources:
        args += ["--source", os.path.join(tmp, f"source_{s}.tsv.gz")]
    args += ["--gold-db", os.path.join(tmp, "gold.sqlite3"), "--alpha", ALPHA,
             "--weight-iterations", "30", "--seed", "7", "--outdir", out,
             "--versions", os.path.join(out, "versions.yml"), "--process-name", "TEST", *extra]
    r = subprocess.run(args, env=ENV, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    return out


def read_classes(path):
    with open(path) as fh:
        return {(r["domain_1"], r["domain_2"]): r["class"] for r in csv.DictReader(fh, delimiter="\t")}


def read_scores(out):
    with gzip.open(os.path.join(out, "ppidm_scores.tsv.gz"), "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_cooccurrence_ignores_same_protein():
    # proteins: 0 {0,1}, 1 {2}, 2 {0}, 3 {0,2}, 4 {3}
    ptr = np.array([0, 2, 3, 4, 6, 7])
    pdom = np.array([0, 1, 2, 0, 0, 2, 3])
    pa, pb = [0, 2, 0, 0], [1, 3, 3, 4]
    for step in (400_000, 1):
        N, n, pc, K = cp.cooccurrence(pa, pb, ptr, pdom, 4, step=step)
        got = {(int(c) // 4, int(c) % 4): int(k) for c, k in zip(pc, K)}
        # (0,1) only from PPI 0-3 across proteins, never from protein 0 alone;
        # (0,2) once per PPI although PPI 0-3 offers it twice; homo (0,0) never
        assert got == {(0, 1): 1, (0, 2): 3, (1, 2): 2, (0, 3): 1, (1, 3): 1}, got
        assert N == 4 and n.tolist() == [4, 3, 3, 1]


def brute_powered(N, n, alpha_low):
    return {(x, y) for x in range(len(n)) for y in range(x + 1, len(n))
            if n[x] and n[y] and hypergeom.cdf(0, N, n[x], n[y]) <= alpha_low}


def test_power_enumeration_matches_brute_force_random():
    rng = np.random.default_rng(0)
    for _ in range(20):
        N = int(rng.integers(20, 300))
        n = rng.integers(0, N // 2, size=30)
        alpha_low = float(10 ** rng.uniform(-8, -1))
        x, y = cp.powered_pairs(N, n, math.log(alpha_low))
        assert set(zip(x.tolist(), y.tolist())) == brute_powered(N, n, alpha_low)
        code = x * len(n) + y
        assert (np.diff(code) > 0).all()  # sorted: lookup() binary-searches it


def test_classes_and_invariants():
    with tempfile.TemporaryDirectory() as tmp:
        write_inputs(tmp)
        out = run(tmp, os.path.join(tmp, "a"))
        pos = read_classes(os.path.join(out, "ppidm_computed.tsv"))
        neg = read_classes(os.path.join(out, "inverse_ppidm_computed.tsv"))
        rep = json.load(open(os.path.join(out, "ppidm_report.json")))

        assert rep["S"] == 5 and rep["h"] == 3
        assert pos[("PF00010", "PF00020")] == "Gold"
        assert pos[("PF00600", "PF00601")] == "Silver"
        assert neg[("PF00100", "PF00200")] == "Gold"
        for pair in (("PF00900", "PF00901"), ("PF00070", "PF00080")):
            assert pair not in pos and pair not in neg
        assert not set(pos) & set(neg)
        assert all(a < b for a, b in list(pos) + list(neg))  # canonical and no homo
        assert rep["gold"]["hetero_in_universe"] == 2
        # gold pairs are emitted and counted, not dropped here
        assert pos[("PF00030", "PF00040")] and rep["classes"]["PPIDM"]["Gold"]["in_3did"] >= 1

        # the same toy data: enumeration == brute force for every source
        prot = toy_proteins()
        acc_idx, ptr, pdom, dom_names = cp.load_proteins(os.path.join(tmp, "protein_pfam.tsv.gz"))
        NP, D = len(acc_idx), len(dom_names)
        assert len(prot) == NP
        rare = (int(np.searchsorted(dom_names, "PF00900")), int(np.searchsorted(dom_names, "PF00901")))
        dense = (int(np.searchsorted(dom_names, "PF00100")), int(np.searchsorted(dom_names, "PF00200")))
        n_dense = 0
        for s in SOURCES:
            keys, _ = cp.load_source(os.path.join(tmp, f"source_{s}.tsv.gz"), acc_idx, NP)
            N, n, _, _ = cp.cooccurrence(keys // NP, keys % NP, ptr, pdom, D)
            x, y = cp.powered_pairs(N, n, math.log(rep["alpha_low"]))
            powered = set(zip(x.tolist(), y.tolist()))
            assert powered == brute_powered(N, n, rep["alpha_low"])
            assert len(x) == rep["sources"][s]["powered_pairs"]
            assert rare not in powered  # absent above for lack of power, not for lack of depletion
            n_dense += dense in powered
        assert n_dense == 4


def test_seed_determinism():
    with tempfile.TemporaryDirectory() as tmp:
        write_inputs(tmp)
        a, b = run(tmp, os.path.join(tmp, "a")), run(tmp, os.path.join(tmp, "b"))
        for f in ("ppidm_computed.tsv", "inverse_ppidm_computed.tsv", "ppidm_scores.tsv.gz", "ppidm_report.json"):
            assert open(os.path.join(a, f), "rb").read() == open(os.path.join(b, f), "rb").read(), f


def test_source_order_does_not_matter():
    """Nextflow stages --source in no guaranteed order; the weight draws must not follow it."""
    with tempfile.TemporaryDirectory() as tmp:
        write_inputs(tmp)
        a = run(tmp, os.path.join(tmp, "a"))
        b = run(tmp, os.path.join(tmp, "b"), sources=list(reversed(SOURCES)))
        for f in ("ppidm_computed.tsv", "inverse_ppidm_computed.tsv", "ppidm_scores.tsv.gz", "ppidm_report.json"):
            assert open(os.path.join(a, f), "rb").read() == open(os.path.join(b, f), "rb").read(), f


def test_pinned_weights_and_threshold():
    w = {"s1": 0.9, "s2": 0.1, "s3": 0.5, "s4": 0.3, "s5": 0.7}
    with tempfile.TemporaryDirectory() as tmp:
        write_inputs(tmp)
        free = read_classes(os.path.join(run(tmp, os.path.join(tmp, "free")), "ppidm_computed.tsv"))
        out = run(tmp, os.path.join(tmp, "pin"), "--weights", ",".join(f"{k}={v}" for k, v in w.items()),
                  "--threshold", "0.5")
        rep = json.load(open(os.path.join(out, "ppidm_report.json")))
        assert rep["weights"] == w and rep["weights_how"] == "pinned"
        assert rep["threshold"] == 0.5 and rep["threshold_how"] == "pinned"
        rows = read_scores(out)
        for r in rows:
            num = sum(w[s] * float(r[f"{s}_C"]) for s in SOURCES)
            den = sum(w[s] for s in SOURCES if float(r[f"{s}_C"]) > 0)
            assert math.isclose(float(r["CS"]), num / den if den else 0.0, rel_tol=1e-5)
            assert (float(r["CS"]) > 0.5) == (r["label"] == "PPIDM")
        pinned = read_classes(os.path.join(out, "ppidm_computed.tsv"))
        assert ("PF00010", "PF00020") in pinned and len(pinned) < len(free)


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"ok {name}")


def test_string_channels_are_one_depletion_vote():
    # The dense never-together pair is depleted in exactly three sources. Under
    # plain names that is 3 >= h = 3 votes (Gold); when two of them are the two
    # STRING channels it is 2 votes (Silver), because they are one network.
    with tempfile.TemporaryDirectory() as tmp:
        write_inputs(tmp)
        edges = toy_edges()
        layouts = {"plain": ["s1", "s2", "s3", "t4", "t5"],
                   "string": ["string_exp", "string_rest", "s3", "t4", "t5"]}
        rows = {0: edges["s1"], 1: edges["s2"], 2: edges["s3"], 3: [], 4: []}
        got = {}
        for layout, names in layouts.items():
            paths = []
            for i, s in enumerate(names):
                path = os.path.join(tmp, layout, f"source_{s}.tsv.gz")
                os.makedirs(os.path.dirname(path), exist_ok=True)
                with gzip.open(path, "wt") as fh:
                    fh.write("acc_a\tacc_b\thi_conf\n")
                    for a, b, hi in rows[i]:
                        a, b = sorted((a, b))
                        fh.write(f"{a}\t{b}\t{hi}\n")
                paths += ["--source", path]
            out = os.path.join(tmp, layout, "out")
            os.makedirs(out)
            r = subprocess.run(
                [sys.executable, SCRIPT, "--protein-pfam", os.path.join(tmp, "protein_pfam.tsv.gz"), *paths,
                 "--gold-db", os.path.join(tmp, "gold.sqlite3"), "--alpha", ALPHA, "--weight-iterations", "5",
                 "--seed", "7", "--outdir", out, "--versions", os.path.join(out, "v.yml"), "--process-name", "T"],
                env=ENV, capture_output=True, text=True)
            assert r.returncode == 0, r.stdout + r.stderr
            got[layout] = read_classes(os.path.join(out, "inverse_ppidm_computed.tsv")).get(("PF00100", "PF00200"))
        assert got == {"plain": "Gold", "string": "Silver"}, got
