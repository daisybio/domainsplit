#!/usr/bin/env python3
"""
Select negative DDIs from a candidate pool (``neg_pool.npz``) for one random
seed, score the selection, and emit the chosen Pfam pairs + a score JSON.

Degree-aware node sampling (DANS, Cappelletti et al. 2024, Bioinformatics
Advances vbae036) matched the negative degree distribution to the positives by
sampling edges with probability proportional to preferential attachment
(PA = deg(a)*deg(b)).  Applied naively to a fixed candidate pool that already has
the positive mean PA, that overshoots: sampling proportional to PA draws edges
with mean PA = E[PA^2]/E[PA] >> pool mean, concentrating on a few hub domains.

This selector keeps the PA-proportional draw but adds a hard per-domain CAP at
the positive degree, and refills in multiple passes until the target negative
count (== positive count) is reached:

  * pass 1 draws ``target`` candidates PA-weighted without replacement and adds
    them while no endpoint exceeds its cap;
  * each refill pass prunes already-picked candidates and any candidate touching
    a saturated domain (it would only be skipped again), renormalises the PA
    weights over the survivors, and draws the deficit (floored at ``min_batch``);
  * the loop stops when ``target`` is hit, or early if the pruned pool can no
    longer supply an eligible edge -- that shortfall is the true feasibility
    ceiling of the fixed pool.

Capping pins the negative degree sequence to the positive one (degree
distribution matched, hub-driven mean-PA blow-up impossible); the PA-weighted
draw fills high-degree domains toward their cap first, reproducing the positive
PA distribution; the refill drives the count to the positive total.

The selection is scored against the positives with a combined objective (lower is
better):  J = w_pa*pa + w_deg*deg + w_cov*cov, where
  pa  = Wasserstein-1 between log1p(PA_neg) and log1p(PA_pos), normalised by the
        spread of the positive log1p(PA);
  deg = Kolmogorov-Smirnov statistic between the per-domain negative degree
        distribution (0 for unused domains) and the positive degree distribution;
  cov = 1 - domains_used_neg / domains_pos.
"""

import argparse
import json

import numpy as np


TAG = "[neg_select]"


def log(msg):
    print(f"{TAG} {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--pool", required=True, help="candidate-pool .npz from BUILD step")
    p.add_argument("--seed", type=int, required=True)
    p.add_argument("--score-out", required=True, help="output score JSON path")
    p.add_argument("--pairs-out", required=True, help="output selected-pairs TSV path")
    p.add_argument("--w-pa", type=float, default=0.5)
    p.add_argument("--w-deg", type=float, default=0.3)
    p.add_argument("--w-cov", type=float, default=0.2)
    p.add_argument(
        "--min-batch",
        type=int,
        default=20,
        help="Minimum candidates drawn per refill pass so the tail keeps progressing.",
    )
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


def select(cand_ai, cand_bi, pa, cap, target, seed, min_batch):
    """PA-weighted, degree-capped, multi-pass refill selection.

    Returns the array of selected candidate indices.
    """
    rng = np.random.default_rng(seed)
    n_cand = cand_ai.size
    remaining = cap.copy()
    picked = np.zeros(n_cand, dtype=bool)
    selected = []

    def draw_and_add(batch):
        elig = np.flatnonzero(
            (~picked) & (remaining[cand_ai] > 0) & (remaining[cand_bi] > 0)
        )
        if elig.size == 0:
            return 0
        w = pa[elig].astype(float)
        total = w.sum()
        p = (w / total) if total > 0 else None
        k = int(min(batch, elig.size))
        chosen = rng.choice(elig, size=k, replace=False, p=p)
        added = 0
        for ci in chosen:
            a = cand_ai[ci]
            b = cand_bi[ci]
            if remaining[a] > 0 and remaining[b] > 0:
                selected.append(int(ci))
                picked[ci] = True
                remaining[a] -= 1
                remaining[b] -= 1
                added += 1
                if len(selected) >= target:
                    break
        return added

    draw_and_add(target)
    n_passes = 1
    while len(selected) < target:
        missing = target - len(selected)
        added = draw_and_add(max(missing, min_batch))
        n_passes += 1
        if added == 0:
            log(f"pool exhausted after {n_passes} passes; "
                f"selected {len(selected)}/{target}")
            break

    log(f"selection done in {n_passes} passes: {len(selected)}/{target} edges")
    return np.array(selected, dtype=np.int64)


def main():
    args = parse_args()

    data = np.load(args.pool, allow_pickle=True)
    cand_a = data["cand_a"]
    cand_b = data["cand_b"]
    pos_dom = data["pos_dom"]
    pos_deg = data["pos_deg"].astype(np.int64)
    pos_edge_pa = data["pos_edge_pa"].astype(np.int64)
    n_positive = int(data["n_positive"])
    n_positive_domains = int(data["n_positive_domains"])

    log(f"pool: {cand_a.size} candidates, {pos_dom.size} positive domains, "
        f"target = {n_positive}")

    # Map every domain to an integer index over the positive-domain universe.
    domain_index = {d: i for i, d in enumerate(pos_dom)}
    cand_ai = np.fromiter((domain_index[a] for a in cand_a), dtype=np.int64,
                          count=cand_a.size)
    cand_bi = np.fromiter((domain_index[b] for b in cand_b), dtype=np.int64,
                          count=cand_b.size)
    cap = pos_deg.copy()
    pa = pos_deg[cand_ai] * pos_deg[cand_bi]

    sel = select(cand_ai, cand_bi, pa, cap, n_positive, args.seed, args.min_batch)

    # --- selected-set statistics ---
    sel_ai = cand_ai[sel]
    sel_bi = cand_bi[sel]
    neg_pa = pa[sel]
    n_sel = int(sel.size)
    mean_pa_neg = float(neg_pa.mean()) if n_sel else 0.0

    neg_deg = np.zeros(pos_dom.size, dtype=np.int64)
    np.add.at(neg_deg, sel_ai, 1)
    np.add.at(neg_deg, sel_bi, 1)
    n_dom = int(np.count_nonzero(neg_deg))

    # --- objective ---
    pos_logpa = np.log1p(pos_edge_pa.astype(float))
    neg_logpa = np.log1p(neg_pa.astype(float))
    spread = float(pos_logpa.max() - pos_logpa.min())
    pa_term = wasserstein1(neg_logpa, pos_logpa) / spread if spread > 0 else 0.0
    deg_term = ks_statistic(neg_deg, pos_deg)
    cov_term = 1.0 - (n_dom / n_positive_domains) if n_positive_domains else 0.0
    j = args.w_pa * pa_term + args.w_deg * deg_term + args.w_cov * cov_term

    pos_mean_pa = float(pos_edge_pa.mean())

    log(f"set=positive n_sel={n_positive} n_dom={n_positive_domains} "
        f"mean_pa={pos_mean_pa:.1f}")
    log(f"set=negative seed={args.seed} J={j:.4f} pa={pa_term:.4f} "
        f"deg={deg_term:.4f} cov={cov_term:.4f} n_sel={n_sel} n_dom={n_dom} "
        f"mean_pa={mean_pa_neg:.1f}")

    score = {
        "seed": int(args.seed),
        "J": j,
        "pa": pa_term,
        "deg": deg_term,
        "cov": cov_term,
        "n_sel": n_sel,
        "n_dom": n_dom,
        "mean_pa": mean_pa_neg,
        "pos_n_sel": n_positive,
        "pos_n_dom": n_positive_domains,
        "pos_mean_pa": pos_mean_pa,
        "w_pa": args.w_pa,
        "w_deg": args.w_deg,
        "w_cov": args.w_cov,
    }
    with open(args.score_out, "w") as fh:
        json.dump(score, fh)

    with open(args.pairs_out, "w") as fh:
        for ci in sel:
            fh.write(f"{cand_a[ci]}\t{cand_b[ci]}\n")


if __name__ == "__main__":
    main()
