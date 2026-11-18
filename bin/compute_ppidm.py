#!/usr/bin/env python3
"""Compute PPIDM (enriched) and INVERSE_PPIDM (depleted) domain pairs from PPI sources.

Statistics as in PLAN_ppidm_in_pipeline.md section 4 (PPIDomainMiner, Algorithm 1
plus the published classification, and its mirror for depletion):

* per source: N PPIs, n(x) PPIs carrying domain x on either protein, K(x, y) PPIs
  with x on one protein and y on the other (hetero pairs only);
* positives from the hi_conf rows: weighted consensus of cosines, threshold
  calibrated on 3did, hypergeometric enrichment per source (Bonferroni over E*);
* negatives from all rows: hypergeometric depletion per source, tested only where
  the source has the power to see one (Bonferroni over every hetero pair of the
  domain universe), enriched nowhere, consensus at or below the threshold.
"""

import argparse
import json
import math
import os
import sqlite3
import sys

import numpy as np
import pandas as pd
import scipy
from scipy.special import gammaln
from scipy.stats import hypergeom, rankdata

from ddi_db_utils import split_sources

BACKGROUND_MAX = 100_000
N_FOLDS = 5
# Depletion votes per *network*, not per source: string_exp and string_rest are
# two channels of one STRING network, and unfiltered both are depleted for
# almost every powered pair (measured: 2,945 of 2,983 INVERSE Gold pairs), so
# counted apart they turned "depleted in 4 of 8" into "in 2 experimental
# databases". Enrichment (enr) and the positive side stay per source.
DEP_VOTE_GROUP = {"string_exp": "string", "string_rest": "string"}
GOLD_SQL = (
    "SELECT da.pfam_id, db.pfam_id, ddi.source, ddi.negative "
    "FROM domain_domain_interaction ddi "
    "JOIN domain da ON da.id = ddi.domain_id_a "
    "JOIN domain db ON db.id = ddi.domain_id_b"
)


def log(msg):
    print(f"[compute_ppidm] {msg}", flush=True)


def fail(msg):
    sys.exit(f"[compute_ppidm] ERROR: {msg}")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--protein-pfam", required=True)
    p.add_argument("--source", action="append", required=True)
    p.add_argument("--gold-db", required=True)
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--weight-iterations", type=int, default=10000)
    p.add_argument("--weights", help="name=w,... pinning every source's weight")
    p.add_argument("--threshold", type=float)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--outdir", default=".")
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


# ------------------------------------------------------------------ inputs
def load_proteins(path):
    """acc -> index map, CSR protein->domain arrays (ptr, pdom) and sorted domain names."""
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    if df["acc"].duplicated().any():
        fail(f"duplicate accessions in {path}")
    ex = df["pfams"].str.split(";").explode()
    ex = ex[ex != ""]
    # sorted names, so index order == string order and canonical pairs agree
    dom_names = np.unique(ex.values)
    D, NP = len(dom_names), len(df)
    occ = np.unique(ex.index.values.astype(np.int64) * D + np.searchsorted(dom_names, ex.values))
    ptr = np.zeros(NP + 1, np.int64)
    ptr[1:] = np.cumsum(np.bincount(occ // D, minlength=NP))
    acc_idx = dict(zip(df["acc"], range(NP)))
    return acc_idx, ptr, (occ % D).astype(np.int64), dom_names


def source_name(path):
    base = os.path.basename(path)
    if not base.startswith("source_"):
        fail(f"source file must be named source_<name>.tsv[.gz]: {path}")
    return base[len("source_"):].split(".")[0]


def load_source(path, acc_idx, NP):
    """Unique canonical PPI keys (all rows, hi_conf rows) over universe proteins."""
    keys, his, skipped = [np.empty(0, np.int64)], [np.empty(0, np.int64)], 0
    # chunked: STRING-rest is ~17M rows, too many Python strings to hold at once
    for df in pd.read_csv(path, sep="\t", dtype={"acc_a": str, "acc_b": str, "hi_conf": np.int64},
                          chunksize=2_000_000):
        a, b = df["acc_a"].map(acc_idx), df["acc_b"].map(acc_idx)
        ok = (a.notna() & b.notna() & (a != b)).values
        skipped += int((~ok).sum())
        a, b = a.values[ok].astype(np.int64), b.values[ok].astype(np.int64)
        key = np.minimum(a, b) * NP + np.maximum(a, b)
        keys.append(key)
        his.append(key[df["hi_conf"].values[ok] == 1])
    if skipped:
        log(f"  {os.path.basename(path)}: {skipped} rows skipped (accession outside universe or self pair)")
    return np.unique(np.concatenate(keys)), np.unique(np.concatenate(his))


def load_gold(path, dom_idx):
    conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    rows = conn.execute(GOLD_SQL).fetchall()
    conn.close()
    names = {tuple(sorted((a, b))) for a, b, src, neg in rows
             if neg is not None and int(neg) == 0 and "3did" in split_sources(src)}
    D = len(dom_idx)
    known = [(dom_idx[a], dom_idx[b]) for a, b in names if a in dom_idx and b in dom_idx]
    codes = np.array(sorted(x * D + y for x, y in known), dtype=np.int64)
    het = codes[codes // D != codes % D]
    n_het_names = sum(a != b for a, b in names)
    return names, codes, het, n_het_names


# ------------------------------------------------------------------ counting
def cooccurrence(pa, pb, ptr, pdom, D, step=400_000):
    """N, n[D], sorted hetero pair codes x*D+y (x<y) and K per code for PPIs pa[i]-pb[i].

    K counts a PPI once per pair, and only across the two proteins: two domains
    on the same protein are not evidence that they interact.
    """
    pa, pb = np.asarray(pa, np.int64), np.asarray(pb, np.int64)
    N = len(pa)
    n = np.zeros(D, np.int64)
    pc, K = np.empty(0, np.int64), np.empty(0, np.int64)

    def side(p):
        lens = ptr[p + 1] - ptr[p]
        off = np.arange(lens.sum()) - np.repeat(np.cumsum(lens) - lens, lens)
        return np.repeat(np.arange(len(p), dtype=np.int64), lens), pdom[np.repeat(ptr[p], lens) + off]

    # chunks partition the PPIs, so per-chunk de-duplication is global de-duplication
    for s in range(0, N, step):
        ia, da = side(pa[s:s + step])
        ib, db = side(pb[s:s + step])
        n += np.bincount(np.unique(np.concatenate([ia * D + da, ib * D + db])) % D, minlength=D)
        m = pd.DataFrame({"i": ia, "x": da}).merge(pd.DataFrame({"i": ib, "y": db}), on="i")
        x, y, i = m["x"].values, m["y"].values, m["i"].values
        ok = x != y
        pair = np.minimum(x, y)[ok] * D + np.maximum(x, y)[ok]
        cc, kc = np.unique(np.unique(i[ok] * (D * D) + pair) % (D * D), return_counts=True)
        pc, inv = np.unique(np.concatenate([pc, cc]), return_inverse=True)
        K = np.bincount(inv, weights=np.concatenate([K, kc]), minlength=len(pc)).astype(np.int64)
    return N, n, pc, K


def lookup(pc, values, codes):
    """values at codes (0 where a code is not in the sorted pc)."""
    if not len(pc):
        return np.zeros(len(codes), values.dtype)
    pos = np.minimum(np.searchsorted(pc, codes), len(pc) - 1)
    return np.where(pc[pos] == codes, values[pos], 0)


def log_p0(N, nx, ny):
    """log hypergeom.cdf(0, N, nx, ny); -inf when nx + ny > N."""
    nx, ny = np.broadcast_arrays(np.asarray(nx, float), np.asarray(ny, float))
    out = np.full(nx.shape, -np.inf)
    ok = nx + ny <= N
    x, y = nx[ok], ny[ok]
    out[ok] = gammaln(N - x + 1) + gammaln(N - y + 1) - gammaln(N + 1) - gammaln(N - x - y + 1)
    return out


def powered_pairs(N, n, log_alpha_low):
    """Exactly the hetero domain pairs (x<y) whose P(K=0) <= alpha_low in this source.

    P0 falls as n_x or n_y grows, so with domains sorted by n descending the
    partners of domain i are a prefix of j > i; binary search its end for every i.
    """
    dom = np.flatnonzero(n > 0)
    order = dom[np.argsort(-n[dom], kind="stable")]
    m, L = n[order], len(order)
    lo, hi = np.arange(L), np.full(L, L)  # power at lo (or lo == i), none at hi
    act = hi - lo > 1
    while act.any():
        mid = (lo + hi) // 2
        p = log_p0(N, m, m[np.minimum(mid, L - 1)]) <= log_alpha_low
        lo = np.where(act & p, mid, lo)
        hi = np.where(act & ~p, mid, hi)
        act = hi - lo > 1
    cnt = lo - np.arange(L)
    ci = np.repeat(np.arange(L), cnt)
    cj = ci + 1 + np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
    x, y = order[ci], order[cj]
    x, y = np.minimum(x, y), np.maximum(x, y)
    o = np.lexsort((y, x))  # sorted by code x*D+y: callers binary-search it
    return x[o], y[o]


# ------------------------------------------------------------------ calibration
def consensus(C, w):
    return (C * w).sum(1) / np.maximum(((C > 0) * w).sum(1), 1e-12)


def auc(pos, neg):
    """Mann-Whitney ROC AUC, average ranks for ties."""
    r = rankdata(np.concatenate([pos, neg]))
    return (r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg))


def prf(pos, neg, ts):
    """F1, precision, recall per threshold t, predicting positive iff score > t."""
    ps, ns = np.sort(pos), np.sort(neg)
    tp = len(ps) - np.searchsorted(ps, ts, side="right")
    fp = len(ns) - np.searchsorted(ns, ts, side="right")
    fn = len(ps) - tp
    f1 = np.where(2 * tp + fn + fp > 0, 2 * tp / np.maximum(2 * tp + fn + fp, 1), 0.0)
    prec = np.where(tp + fp > 0, tp / np.maximum(tp + fp, 1), 0.0)
    rec = tp / max(len(ps), 1)
    return f1, prec, rec


def grid(lo, hi, step):
    k0, k1 = math.ceil(max(lo, 0.0) / step - 1e-9), math.floor(min(hi, 1.0) / step + 1e-9)
    return np.round(np.arange(k0, k1 + 1) * step, 6)


def main():
    args = parse_args()
    out = lambda name: os.path.join(args.outdir, name)  # noqa: E731
    # weight draws map to sources by position, so fix the order independent of staging
    args.source = sorted(args.source, key=source_name)
    names = [source_name(p) for p in args.source]
    if len(set(names)) != len(names):
        fail(f"duplicate source names: {names}")
    S, h = len(names), math.ceil(len(names) / 2)

    acc_idx, ptr, pdom, dom_names = load_proteins(args.protein_pfam)
    NP, D = len(acc_idx), len(dom_names)
    dom_idx = {d: i for i, d in enumerate(dom_names)}
    log(f"universe: {NP} proteins, {D} domains; {S} sources, h = {h}")

    POS, NEG = {}, {}
    for name, path in zip(names, args.source):
        keys_all, keys_hi = load_source(path, acc_idx, NP)
        NEG[name] = cooccurrence(keys_all // NP, keys_all % NP, ptr, pdom, D)
        POS[name] = cooccurrence(keys_hi // NP, keys_hi % NP, ptr, pdom, D)
        log(f"  {name}: N all {len(keys_all)} (pairs {len(NEG[name][2])}), "
            f"N hi_conf {len(keys_hi)} (pairs {len(POS[name][2])})")

    gold_names, gold_all, gold_het, n_het_names = load_gold(args.gold_db, dom_idx)
    log(f"3did gold: {len(gold_names)} pairs, {n_het_names} hetero, "
        f"{len(gold_het)} hetero with both domains in the universe")
    if not len(gold_het):
        fail("no hetero 3did pair with both domains in the universe in --gold-db")

    # ---------------------------------------------------------- consensus
    allc = np.unique(np.concatenate([POS[s][2] for s in names]))
    if not len(allc):
        fail("no co-occurring domain pair in any hi_conf source")
    a, b = allc // D, allc % D
    C = np.zeros((len(allc), S))
    Kpos = {}
    for j, s in enumerate(names):
        N, n, pc, K = POS[s]
        Kpos[s] = lookup(pc, K, allc)
        C[:, j] = Kpos[s] / np.sqrt(np.maximum(n[a] * n[b], 1))

    isg = np.isin(allc, gold_het)
    gidx, n_unscored_gold = np.flatnonzero(isg), len(gold_het) - int(isg.sum())
    if args.weights:
        pinned = dict(kv.split("=") for kv in args.weights.split(","))
        if set(pinned) != set(names):
            fail(f"--weights must pin exactly the sources {names}, got {sorted(pinned)}")
        w_best = np.array([float(pinned[s]) for s in names])
        if (w_best <= 0).any():
            fail("--weights must be > 0")
        auc_best = auc_eq = None
        log(f"weights pinned: {dict(zip(names, w_best.tolist()))}")
    else:
        rng_w = np.random.default_rng(args.seed)
        pool = np.flatnonzero(~isg)
        if not len(pool):
            fail("no scored non-gold pair to serve as AUC background")
        bgi = np.sort(rng_w.choice(pool, min(BACKGROUND_MAX, len(pool)), replace=False))
        Cg, Cb, zeros = C[gidx], C[bgi], np.zeros(n_unscored_gold)
        score = lambda w: auc(np.concatenate([consensus(Cg, w), zeros]), consensus(Cb, w))  # noqa: E731
        w_best = np.ones(S)
        auc_best = auc_eq = score(w_best)
        for _ in range(args.weight_iterations):
            w = rng_w.integers(1, 101, S) / 100
            v = score(w)
            if v > auc_best:
                w_best, auc_best = w, v
        log(f"weights (AUC {auc_best:.4f}, equal weights {auc_eq:.4f}, background {len(bgi)}): "
            + ", ".join(f"{s}={w:.2f}" for s, w in zip(names, w_best)))
    CS = consensus(C, w_best)

    # ---------------------------------------------------------- threshold
    rng_t = np.random.default_rng(args.seed)
    gs = lookup(allc, CS, gold_het)
    ga, gb = gold_het // D, gold_het % D
    sb = rng_t.permutation(gb)
    sc = np.minimum(ga, sb) * D + np.maximum(ga, sb)
    sc = np.unique(sc[(ga != sb) & ~np.isin(sc, gold_all)])
    ns = lookup(allc, CS, sc)

    def split(x):
        perm = rng_t.permutation(len(x))
        return x[perm[len(x) // 3:]], x[perm[:len(x) // 3]]

    pos_tr, pos_te = split(gs)
    neg_tr, neg_te = split(ns)
    fold_p = rng_t.permutation(len(pos_tr)) % N_FOLDS
    fold_n = rng_t.permutation(len(neg_tr)) % N_FOLDS

    def cv_f1(ts):
        return np.mean([prf(pos_tr[fold_p == k], neg_tr[fold_n == k], ts)[0] for k in range(N_FOLDS)], axis=0)

    if args.threshold is not None:
        T, how = float(args.threshold), "pinned"
    else:
        ts = grid(0.0, 1.0, 0.01)
        T = float(ts[np.argmax(cv_f1(ts))])  # argmax takes the first, i.e. smallest, on ties
        for prev, step in ((0.01, 0.001), (0.001, 0.0001)):
            ts = grid(T - prev, T + prev, step)
            T = float(ts[np.argmax(cv_f1(ts))])
        how = "scan"
    f1_cv = float(cv_f1(np.array([T]))[0])
    f1_te, p_te, r_te = (float(v[0]) for v in prf(pos_te, neg_te, np.array([T])))
    log(f"threshold {T:.4f} ({how}): CV F1 {f1_cv:.3f}; test F1 {f1_te:.3f} precision {p_te:.3f} recall {r_te:.3f} "
        f"(shuffle negatives {len(sc)}, positives {len(gold_het)})")

    # ---------------------------------------------------------- positives
    E = CS > T
    alpha_up = args.alpha / max(int(E.sum()), 1)
    ec, sig = allc[E], np.zeros(int(E.sum()), int)
    for s in names:
        N, n, _, _ = POS[s]
        k = Kpos[s][E]
        has = k > 0
        pu = np.ones(len(k))
        pu[has] = hypergeom.sf(k[has] - 1, N, n[ec[has] // D], n[ec[has] % D])
        sig += has & (pu <= alpha_up)
    pos_class = np.where(sig >= h, "Gold", np.where(sig > 0, "Silver", "Bronze"))

    # ---------------------------------------------------------- negatives
    present = np.zeros(D, bool)
    for s in names:
        present |= NEG[s][1] > 0
    Du = int(present.sum())
    alpha_low = args.alpha / (Du * (Du - 1) / 2) if Du > 1 else 0.0
    log_alpha_low = math.log(alpha_low) if alpha_low > 0 else -np.inf
    plow, src_stats = {}, {}
    for s in names:
        N, n, pc, K = NEG[s]
        x, y = powered_pairs(N, n, log_alpha_low)
        code = x * D + y
        pl = hypergeom.cdf(lookup(pc, K, code), N, n[x], n[y]) if len(code) else np.empty(0)
        plow[s] = (code, pl)
        src_stats[s] = {"N_hi_conf": POS[s][0], "N_all": N,
                        "domains_hi_conf": int((POS[s][1] > 0).sum()), "domains_all": int((n > 0).sum()),
                        "cooccurring_pairs_hi_conf": len(POS[s][2]), "cooccurring_pairs_all": len(pc),
                        "powered_pairs": len(code), "depleted_pairs": int((pl <= alpha_low).sum())}
        log(f"  {s}: powered {len(code)}, depleted {src_stats[s]['depleted_pairs']}")
    cand = np.unique(np.concatenate([plow[s][0] for s in names]))
    enr = np.zeros(len(cand), int)
    dep_by_group = {}
    for s in names:
        code, pl = plow[s]
        g = DEP_VOTE_GROUP.get(s, s)
        dep_by_group[g] = dep_by_group.get(g, np.zeros(len(cand), bool)) | (lookup(code, (pl <= alpha_low).astype(int), cand) > 0)
        N, n, pc, K = NEG[s]
        k = lookup(pc, K, cand)
        has = k > 0
        pu = np.ones(len(cand))
        pu[has] = hypergeom.sf(k[has] - 1, N, n[cand[has] // D], n[cand[has] % D])
        enr += has & (pu <= alpha_up)
    dep = sum(v.astype(int) for v in dep_by_group.values()) if dep_by_group else np.zeros(len(cand), int)
    ok = (enr == 0) & (lookup(allc, CS, cand) <= T)
    inv_codes = cand[ok & (dep > 0)]
    inv_class = np.where(dep[ok & (dep > 0)] >= h, "Gold", "Silver")
    assert not np.intersect1d(ec, inv_codes).size, "PPIDM and INVERSE_PPIDM overlap"
    assert not ((ec // D == ec % D).any() or (inv_codes // D == inv_codes % D).any()), "homo pair emitted"

    # ---------------------------------------------------------- outputs
    rank = {"Gold": 0, "Silver": 1, "Bronze": 2}
    frames = []
    for label, codes, cls, fname in (("PPIDM", ec, pos_class, "ppidm_computed.tsv"),
                                     ("INVERSE_PPIDM", inv_codes, inv_class, "inverse_ppidm_computed.tsv")):
        df = pd.DataFrame({"code": codes, "domain_1": dom_names[codes // D], "domain_2": dom_names[codes % D],
                           "label": label, "class": cls, "rank": [rank[c] for c in cls]})
        df = df.sort_values(["rank", "domain_1", "domain_2"], kind="stable")
        df[["domain_1", "domain_2", "class"]].to_csv(out(fname), sep="\t", index=False)
        frames.append(df)
    em = pd.concat(frames, ignore_index=True)
    codes = em["code"].values
    em["CS"] = lookup(allc, CS, codes)
    for s in names:
        N, n, pc, K = POS[s]
        k = lookup(pc, K, codes)
        nx, ny = n[codes // D], n[codes % D]
        has = k > 0
        pu = np.full(len(codes), np.nan)
        pu[has] = hypergeom.sf(k[has] - 1, N, nx[has], ny[has])
        pcode, pl = plow[s]
        pos = np.minimum(np.searchsorted(pcode, codes), max(len(pcode) - 1, 0))
        powered = (pcode[pos] == codes) if len(pcode) else np.zeros(len(codes), bool)
        em[f"{s}_K"], em[f"{s}_nx"], em[f"{s}_ny"] = k, nx, ny
        em[f"{s}_C"] = k / np.sqrt(np.maximum(nx * ny, 1))
        em[f"{s}_pup"] = pu
        em[f"{s}_plow"] = np.where(powered, pl[pos] if len(pl) else np.nan, np.nan)
    em.drop(columns=["code", "rank"]).to_csv(
        out("ppidm_scores.tsv.gz"), sep="\t", index=False, na_rep="NA", float_format="%.6g",
        compression={"method": "gzip", "mtime": 0})  # mtime 0: byte-identical reruns

    def counts(codes, cls, classes):
        return {c: {"pairs": int((cls == c).sum()), "in_3did": int(np.isin(codes[cls == c], gold_all).sum())}
                for c in classes}

    report = {
        "S": S, "h": h, "seed": args.seed, "alpha": args.alpha,
        "depletion_vote_groups": {s: DEP_VOTE_GROUP.get(s, s) for s in names},
        "sources": src_stats,
        "weights": dict(zip(names, w_best.tolist())),
        "weights_how": "pinned" if args.weights else "search",
        "weight_iterations": 0 if args.weights else args.weight_iterations,
        "weight_auc_best": auc_best, "weight_auc_equal": auc_eq,
        "threshold": T, "threshold_how": how,
        "threshold_scan": "0.01 over [0,1], then 0.001 and 0.0001 around the previous best "
                          "(deviation from the paper's fixed phase windows: measured optimum lies outside them)",
        "threshold_cv_f1": f1_cv, "threshold_test_f1": f1_te,
        "threshold_test_precision": p_te, "threshold_test_recall": r_te,
        "shuffle": {"positives": len(gold_het), "negatives": len(sc),
                    "train_positives": len(pos_tr), "train_negatives": len(neg_tr),
                    "test_positives": len(pos_te), "test_negatives": len(neg_te)},
        "gold": {"pairs": len(gold_names), "hetero": n_het_names, "hetero_in_universe": len(gold_het),
                 "hetero_in_universe_unscored": n_unscored_gold},
        "scored_pairs": len(allc), "E_star": int(E.sum()), "alpha_up": alpha_up,
        "domain_universe": Du, "alpha_low": alpha_low, "powered_pairs_any_source": len(cand),
        "classes": {"PPIDM": counts(ec, pos_class, ("Gold", "Silver", "Bronze")),
                    "INVERSE_PPIDM": counts(inv_codes, inv_class, ("Gold", "Silver"))},
    }
    with open(out("ppidm_report.json"), "w") as fh:
        json.dump(report, fh, indent=2)
    log("PPIDM " + ", ".join(f"{c} {v['pairs']} (3did {v['in_3did']})" for c, v in report["classes"]["PPIDM"].items())
        + " | INVERSE_PPIDM "
        + ", ".join(f"{c} {v['pairs']} (3did {v['in_3did']})" for c, v in report["classes"]["INVERSE_PPIDM"].items()))

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    numpy: {np.__version__}\n    scipy: {scipy.__version__}\n    pandas: {pd.__version__}\n")


if __name__ == "__main__":
    main()
