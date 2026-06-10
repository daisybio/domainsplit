#!/usr/bin/env python3
"""Local unit-check for bin/select_ppi_negative_dans.py (no Nextflow, no cluster).

Builds a tiny synthetic candidate pool and runs both DANS methods, asserting the
core invariants:

  * deletion        -- draws exactly n_positive_r edges, all from the candidate
                       pool, all endpoints within the pool-domain universe.
  * random_addition -- draws exactly n_positive edges, none a positive/forbidden
                       pair, no self-pairs, no duplicates, and reaches domains
                       that are absent from the candidate pool.

Run directly (`python3 tests/python/test_select_ppi_negative_dans.py`) or via
pytest.
"""

import os
import subprocess
import sys
import tempfile

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
SELECTOR = os.path.join(BIN, "select_ppi_negative_dans.py")


def pf(i):
    return f"PF{i:05d}"


def build_pool(path):
    # Full 3did positive edges (canonical ascending).
    pos_edges = [(1, 2), (1, 3), (1, 4), (2, 3), (5, 6),
                 (7, 8), (9, 10), (1, 9), (2, 10)]
    # Candidate pool: fresh (non-positive) pairs over domains 1..6.
    cand = [(1, 5), (1, 6), (2, 5), (2, 6), (3, 5),
            (3, 6), (4, 5), (4, 6), (3, 4)]
    pool_domains = sorted({d for e in cand for d in e})        # 1..6

    # Reduced positives = positive edges with both endpoints in the pool.
    pos_r = [(a, b) for a, b in pos_edges
             if a in pool_domains and b in pool_domains]
    deg_r = {d: 0 for d in pool_domains}
    for a, b in pos_r:
        deg_r[a] += 1
        deg_r[b] += 1
    pool_dom = [pf(d) for d in pool_domains]
    pool_deg_r = np.array([deg_r[d] for d in pool_domains], dtype=np.int64)
    pos_edge_pa_r = np.array([deg_r[a] * deg_r[b] for a, b in pos_r], dtype=np.int64)

    # Full positive degrees over all 10 domains.
    all_dom = list(range(1, 11))
    deg = {d: 0 for d in all_dom}
    for a, b in pos_edges:
        deg[a] += 1
        deg[b] += 1
    pos_dom = [pf(d) for d in all_dom]
    pos_deg = np.array([deg[d] for d in all_dom], dtype=np.int64)
    pos_edge_pa = np.array([deg[a] * deg[b] for a, b in pos_edges], dtype=np.int64)

    np.savez(
        path,
        cand_a=np.array([pf(a) for a, b in cand], dtype=object),
        cand_b=np.array([pf(b) for a, b in cand], dtype=object),
        pool_dom=np.array(pool_dom, dtype=object),
        pool_deg_r=pool_deg_r,
        pos_edge_pa_r=pos_edge_pa_r,
        n_positive_r=np.int64(len(pos_r)),
        n_positive_domains_r=np.int64(sum(1 for d in pool_domains if deg_r[d])),
        pos_a=np.array([pf(a) for a, b in pos_edges], dtype=object),
        pos_b=np.array([pf(b) for a, b in pos_edges], dtype=object),
        pos_dom=np.array(pos_dom, dtype=object),
        pos_deg=pos_deg,
        pos_edge_pa=pos_edge_pa,
        n_positive=np.int64(len(pos_edges)),
        n_positive_domains=np.int64(len(all_dom)),
        forbidden_a=np.array([pf(a) for a, b in pos_edges], dtype=object),
        forbidden_b=np.array([pf(b) for a, b in pos_edges], dtype=object),
    )
    return {
        "cand": {tuple(sorted((pf(a), pf(b)))) for a, b in cand},
        "pool_domains": {pf(d) for d in pool_domains},
        "n_positive_r": len(pos_r),
        "n_positive": len(pos_edges),
        "forbidden": {(pf(a), pf(b)) for a, b in pos_edges},
        "pool_only": {pf(d) for d in pool_domains},
        "extra_domains": {pf(9), pf(10)},
    }


def run_method(pool_path, method, workdir):
    score = os.path.join(workdir, f"score_{method}.json")
    pairs = os.path.join(workdir, f"pairs_{method}.tsv")
    env = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))
    subprocess.run(
        [sys.executable, SELECTOR, "--pool", pool_path, "--method", method,
         "--seed", "7", "--score-out", score, "--pairs-out", pairs],
        check=True, env=env,
    )
    out = []
    with open(pairs) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line:
                a, b = line.split("\t")
                out.append((a, b))
    return out


def test_dans_methods():
    with tempfile.TemporaryDirectory() as tmp:
        pool_path = os.path.join(tmp, "neg_pool.npz")
        meta = build_pool(pool_path)

        # --- deletion ---
        del_pairs = run_method(pool_path, "deletion", tmp)
        assert len(del_pairs) == meta["n_positive_r"], \
            f"deletion count {len(del_pairs)} != {meta['n_positive_r']}"
        for a, b in del_pairs:
            assert tuple(sorted((a, b))) in meta["cand"], f"{a},{b} not in pool"
            assert a in meta["pool_domains"] and b in meta["pool_domains"]
        assert len({tuple(sorted(p)) for p in del_pairs}) == len(del_pairs), "dup in deletion"

        # --- random_addition ---
        rand_pairs = run_method(pool_path, "random_addition", tmp)
        assert len(rand_pairs) == meta["n_positive"], \
            f"random_addition count {len(rand_pairs)} != {meta['n_positive']}"
        seen = set()
        used_domains = set()
        for a, b in rand_pairs:
            assert a != b, f"self pair {a}"
            key = (a, b) if a <= b else (b, a)
            assert key not in meta["forbidden"], f"{key} is a positive/forbidden pair"
            assert key not in seen, f"duplicate {key}"
            seen.add(key)
            used_domains.update((a, b))
        # DANS over the full positive set must be able to reach domains outside
        # the candidate pool.
        assert used_domains & meta["extra_domains"], \
            "random_addition never reached the pool-absent domains"

    print("OK: both DANS methods satisfy invariants")


if __name__ == "__main__":
    test_dans_methods()
