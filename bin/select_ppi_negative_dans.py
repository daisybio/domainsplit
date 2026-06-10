#!/usr/bin/env python3
"""
Select negative DDIs by Degree-Aware Node Sampling (DANS, Cappelletti et al.
2024, Bioinformatics Advances vbae036) for one of two methods, and write the
chosen Pfam pairs plus a small score JSON (reporting only -- there is no
multi-seed pick-best).

DANS draws a negative edge by sampling two endpoints proportional to node degree
(equivalently: take the source of one random positive edge and the destination
of another) and accepts the pair iff it is not an existing edge.  It is
UNCAPPED: the negative degree *distribution* tracks the positive one without
pinning an exact per-node degree sequence.

  * method "deletion"        -- DANS restricted to the PPI-derived candidate pool
    (``cand_a``/``cand_b``), with the positive degrees first reduced to the
    candidate-domain universe (``pool_deg_r``).  Candidate edges are drawn
    without replacement with probability proportional to the reduced preferential
    attachment PA_r = deg_r(a)*deg_r(b); target = ``n_positive_r``.

  * method "random_addition" -- plain DANS over the *full* positive set: sample
    node-pairs from the positive endpoint multiset (``pos_a``/``pos_b``),
    rejecting self-pairs, existing edges (``forbidden_*``) and duplicates;
    target = ``n_positive``.  Domains absent from the candidate pool are reachable
    here, so coverage and the degree distribution match the full positives.

The selection is scored against the method-appropriate positives with a combined
objective (lower is better), reported for inspection only:
  J = w_pa*pa + w_deg*deg + w_cov*cov, where
  pa  = Wasserstein-1 between log1p(PA_neg) and log1p(PA_pos), normalised by the
        spread of the positive log1p(PA);
  deg = Kolmogorov-Smirnov statistic between the per-domain negative degree
        distribution (0 for unused domains) and the positive degree distribution;
  cov = 1 - domains_used_neg / domains_pos.
"""

import argparse
import json

import numpy as np

from ddi_db_utils import pfam_sort_key


TAG = "[neg_select]"


def log(msg):
    print(f"{TAG} {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--pool", required=True, help="candidate-pool .npz from BUILD step")
    p.add_argument("--method", required=True,
                   choices=["deletion", "random_addition"])
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--score-out", required=True, help="output score JSON path")
    p.add_argument("--pairs-out", required=True, help="output selected-pairs TSV path")
    p.add_argument("--w-pa", type=float, default=0.5)
    p.add_argument("--w-deg", type=float, default=0.3)
    p.add_argument("--w-cov", type=float, default=0.2)
    return p.parse_args()


def wasserstein1(x, y):
    """1-Wasserstein distance between two 1D empirical samples (numpy only)."""
    x = np.sort(np.asarray(x, dtype=float))
    y = np.sort(np.asarray(y, dtype=float))
    grid = np.concatenate([x, y])
    grid.sort()
    cx = np.searchsorted(x, grid, side="right") / x.size
    cy = np.searchsorted(y, grid, side="right") / y.size
    deltas = np.diff(grid)
    return float(np.sum(np.abs(cx[:-1] - cy[:-1]) * deltas))


def ks_statistic(x, y):
    """Two-sample Kolmogorov-Smirnov statistic (numpy only)."""
    x = np.sort(np.asarray(x, dtype=float))
    y = np.sort(np.asarray(y, dtype=float))
    grid = np.concatenate([x, y])
    grid.sort()
    cx = np.searchsorted(x, grid, side="right") / x.size
    cy = np.searchsorted(y, grid, side="right") / y.size
    return float(np.max(np.abs(cx - cy)))


def select_deletion(cand_a, cand_b, pool_dom, pool_deg_r, target, seed):
    """DANS over the fixed candidate pool: draw `target` edges without
    replacement with probability proportional to the reduced PA. No cap.

    Returns (selected_indices, cand_ai, cand_bi) where cand_a*/cand_b* index
    into pool_dom.
    """
    rng = np.random.default_rng(seed)
    domain_index = {d: i for i, d in enumerate(pool_dom)}
    cand_ai = np.fromiter((domain_index[a] for a in cand_a), dtype=np.int64,
                          count=cand_a.size)
    cand_bi = np.fromiter((domain_index[b] for b in cand_b), dtype=np.int64,
                          count=cand_b.size)
    pa = (pool_deg_r[cand_ai] * pool_deg_r[cand_bi]).astype(float)
    n_cand = int(cand_a.size)
    k = int(min(target, n_cand))
    if k < target:
        log(f"deletion: candidate pool smaller than target "
            f"({n_cand} < {target}); taking all")
    total = pa.sum()
    p = (pa / total) if total > 0 else None
    sel = (rng.choice(n_cand, size=k, replace=False, p=p)
           if k > 0 else np.empty(0, dtype=np.int64))
    log(f"deletion: selected {sel.size}/{target} candidate edges")
    return sel, cand_ai, cand_bi


def select_random_addition(pos_a, pos_b, forbidden, target, seed):
    """Canonical DANS over the full positive set: sample node-pairs from the
    positive endpoint multiset (so endpoints are drawn proportional to degree),
    rejecting self-pairs, existing edges and duplicates. No cap.
    """
    rng = np.random.default_rng(seed)
    endpoints = np.concatenate([pos_a, pos_b])
    m = int(endpoints.size)
    picked = set()
    out_a = []
    out_b = []
    attempts = 0
    max_attempts = 200 * target + 1000
    while len(out_a) < target and attempts < max_attempts:
        need = target - len(out_a)
        batch = int(min(max(need * 2, 1024), 5_000_000))
        ui = rng.integers(0, m, size=batch)
        vi = rng.integers(0, m, size=batch)
        attempts += batch
        for iu, iv in zip(ui.tolist(), vi.tolist()):
            a = endpoints[iu]
            b = endpoints[iv]
            if a == b:
                continue
            key = (a, b) if pfam_sort_key(a) <= pfam_sort_key(b) else (b, a)
            if key in forbidden or key in picked:
                continue
            picked.add(key)
            out_a.append(key[0])
            out_b.append(key[1])
            if len(out_a) >= target:
                break
    if len(out_a) < target:
        log(f"random_addition: only {len(out_a)}/{target} edges after "
            f"{attempts} attempts (forbidden/duplicate saturation)")
    else:
        log(f"random_addition: selected {target}/{target} edges in "
            f"{attempts} attempts")
    return np.array(out_a, dtype=object), np.array(out_b, dtype=object)


def main():
    args = parse_args()
    data = np.load(args.pool, allow_pickle=True)

    if args.method == "deletion":
        cand_a = data["cand_a"]
        cand_b = data["cand_b"]
        pool_dom = data["pool_dom"]
        dom_deg = data["pool_deg_r"].astype(np.int64)
        pos_edge_pa = data["pos_edge_pa_r"].astype(np.int64)
        target = int(data["n_positive_r"])
        n_pos_domains = int(data["n_positive_domains_r"])

        log(f"pool: {cand_a.size} candidate edges, {pool_dom.size} pool domains, "
            f"target = {target}")
        sel, cand_ai, cand_bi = select_deletion(
            cand_a, cand_b, pool_dom, dom_deg, target, args.seed
        )
        neg_ai = cand_ai[sel]
        neg_bi = cand_bi[sel]
        out_a = cand_a[sel]
        out_b = cand_b[sel]
    else:
        pos_a = data["pos_a"]
        pos_b = data["pos_b"]
        pos_dom = data["pos_dom"]
        dom_deg = data["pos_deg"].astype(np.int64)
        pos_edge_pa = data["pos_edge_pa"].astype(np.int64)
        target = int(data["n_positive"])
        n_pos_domains = int(data["n_positive_domains"])
        forbidden = set(zip(data["forbidden_a"].tolist(),
                            data["forbidden_b"].tolist()))

        log(f"positives: {pos_a.size} edges, {pos_dom.size} domains, "
            f"{len(forbidden)} forbidden pairs, target = {target}")
        out_a, out_b = select_random_addition(
            pos_a, pos_b, forbidden, target, args.seed
        )
        domain_index = {d: i for i, d in enumerate(pos_dom)}
        neg_ai = np.fromiter((domain_index[a] for a in out_a), dtype=np.int64,
                             count=out_a.size)
        neg_bi = np.fromiter((domain_index[b] for b in out_b), dtype=np.int64,
                             count=out_b.size)

    # --- selected-set statistics / objective (reporting only) ---
    n_sel = int(out_a.size)
    neg_pa = (dom_deg[neg_ai] * dom_deg[neg_bi]).astype(np.int64)
    mean_pa_neg = float(neg_pa.mean()) if n_sel else 0.0

    neg_deg = np.zeros(dom_deg.size, dtype=np.int64)
    np.add.at(neg_deg, neg_ai, 1)
    np.add.at(neg_deg, neg_bi, 1)
    n_dom = int(np.count_nonzero(neg_deg))

    pos_logpa = np.log1p(pos_edge_pa.astype(float))
    neg_logpa = np.log1p(neg_pa.astype(float))
    spread = float(pos_logpa.max() - pos_logpa.min()) if pos_logpa.size else 0.0
    pa_term = wasserstein1(neg_logpa, pos_logpa) / spread if spread > 0 else 0.0
    deg_term = ks_statistic(neg_deg, dom_deg)
    cov_term = 1.0 - (n_dom / n_pos_domains) if n_pos_domains else 0.0
    j = args.w_pa * pa_term + args.w_deg * deg_term + args.w_cov * cov_term
    pos_mean_pa = float(pos_edge_pa.mean()) if pos_edge_pa.size else 0.0

    log(f"method={args.method} set=positive n_sel={target} n_dom={n_pos_domains} "
        f"mean_pa={pos_mean_pa:.1f}")
    log(f"method={args.method} set=negative seed={args.seed} J={j:.4f} "
        f"pa={pa_term:.4f} deg={deg_term:.4f} cov={cov_term:.4f} n_sel={n_sel} "
        f"n_dom={n_dom} mean_pa={mean_pa_neg:.1f}")

    score = {
        "method": args.method,
        "seed": int(args.seed),
        "J": j,
        "pa": pa_term,
        "deg": deg_term,
        "cov": cov_term,
        "n_sel": n_sel,
        "n_dom": n_dom,
        "mean_pa": mean_pa_neg,
        "pos_n_sel": target,
        "pos_n_dom": n_pos_domains,
        "pos_mean_pa": pos_mean_pa,
        "w_pa": args.w_pa,
        "w_deg": args.w_deg,
        "w_cov": args.w_cov,
    }
    with open(args.score_out, "w") as fh:
        json.dump(score, fh)

    with open(args.pairs_out, "w") as fh:
        for a, b in zip(out_a.tolist(), out_b.tolist()):
            fh.write(f"{a}\t{b}\n")


if __name__ == "__main__":
    main()
