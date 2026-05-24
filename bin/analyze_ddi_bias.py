#!/usr/bin/env python3
"""
Analyze domain-domain interaction dataset for biases that can inflate
downstream prediction performance (Bernett et al. 2024).

Checks:
  1. Degree distribution mismatch between positive and negative DDIs
  2. Degree ratio analysis (fraction of interactions that are positive)
  3. Hub / study bias (overrepresented domains)

Outputs JSON report + matplotlib plots to --outdir.
"""

import argparse
import json
import sqlite3
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--db", required=True, help="Path to enriched domainsplit.sqlite3")
    p.add_argument("--outdir", required=True, help="Output directory for reports and plots")
    return p.parse_args()


def compute_degree(conn, negative_flag):
    """Compute per-domain degree for positive (0) or negative (1) DDIs."""
    rows = conn.execute(
        "SELECT domain_id_a, domain_id_b FROM domain_domain_interaction WHERE negative = ?",
        (negative_flag,),
    ).fetchall()
    degree = Counter()
    for a, b in rows:
        degree[a] += 1
        degree[b] += 1
    return degree


def ks_test(sample_a, sample_b):
    """Two-sample Kolmogorov-Smirnov test (no scipy dependency)."""
    a = np.sort(sample_a)
    b = np.sort(sample_b)
    all_vals = np.sort(np.concatenate([a, b]))
    cdf_a = np.searchsorted(a, all_vals, side="right") / len(a)
    cdf_b = np.searchsorted(b, all_vals, side="right") / len(b)
    ks_stat = float(np.max(np.abs(cdf_a - cdf_b)))
    n = len(a) * len(b) / (len(a) + len(b))
    # Approximate p-value using asymptotic formula
    lam = (np.sqrt(n) + 0.12 + 0.11 / np.sqrt(n)) * ks_stat
    if lam < 1e-10:
        p_value = 1.0
    else:
        # Kolmogorov distribution approximation
        p_value = 2.0 * sum(
            ((-1) ** (k - 1)) * np.exp(-2.0 * k * k * lam * lam)
            for k in range(1, 101)
        )
        p_value = max(0.0, min(1.0, p_value))
    return ks_stat, p_value


def gini_coefficient(values):
    """Compute Gini coefficient of a distribution."""
    arr = np.sort(np.array(values, dtype=float))
    n = len(arr)
    if n == 0 or arr.sum() == 0:
        return 0.0
    index = np.arange(1, n + 1)
    return float((2.0 * np.sum(index * arr) - (n + 1) * np.sum(arr)) / (n * np.sum(arr)))


def plot_degree_distributions(pos_degrees, neg_degrees, outdir):
    """Overlay histogram of positive vs negative degree distributions."""
    fig, ax = plt.subplots(figsize=(8, 5))
    all_vals = list(pos_degrees.values()) + list(neg_degrees.values())
    if not all_vals:
        return
    max_deg = max(all_vals)
    bins = np.logspace(0, np.log10(max_deg + 1), 30)

    ax.hist(list(pos_degrees.values()), bins=bins, alpha=0.6, label="Positive", color="steelblue")
    ax.hist(list(neg_degrees.values()), bins=bins, alpha=0.6, label="Negative", color="coral")
    ax.set_xscale("log")
    ax.set_xlabel("Domain degree")
    ax.set_ylabel("Count")
    ax.set_title("Degree distribution: positive vs negative DDIs")
    ax.legend()
    fig.tight_layout()
    fig.savefig(Path(outdir) / "degree_distribution.png", dpi=150)
    plt.close(fig)


def plot_degree_ratios(degree_ratios, outdir):
    """Histogram of per-domain degree ratios."""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(degree_ratios, bins=20, range=(0, 1), color="mediumpurple", edgecolor="black", alpha=0.8)
    ax.set_xlabel("Degree ratio (n_pos / n_total)")
    ax.set_ylabel("Number of domains")
    ax.set_title("Per-domain degree ratio distribution")
    ax.axvline(0.5, color="red", linestyle="--", alpha=0.5, label="Balanced (0.5)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(Path(outdir) / "degree_ratios.png", dpi=150)
    plt.close(fig)


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    conn = sqlite3.connect(args.db)
    report = {}

    # --- 1. Degree distributions ---
    pos_deg = compute_degree(conn, 0)
    neg_deg = compute_degree(conn, 1)

    pos_vals = np.array(list(pos_deg.values()), dtype=float) if pos_deg else np.array([0.0])
    neg_vals = np.array(list(neg_deg.values()), dtype=float) if neg_deg else np.array([0.0])

    ks_stat, ks_p = ks_test(pos_vals, neg_vals) if len(pos_vals) > 1 and len(neg_vals) > 1 else (0.0, 1.0)

    report["degree_distribution"] = {
        "n_domains_in_positive": len(pos_deg),
        "n_domains_in_negative": len(neg_deg),
        "pos_degree_mean": float(np.mean(pos_vals)),
        "pos_degree_median": float(np.median(pos_vals)),
        "neg_degree_mean": float(np.mean(neg_vals)),
        "neg_degree_median": float(np.median(neg_vals)),
        "ks_statistic": ks_stat,
        "ks_p_value": ks_p,
        "interpretation": "p < 0.05 indicates significantly different distributions"
    }

    plot_degree_distributions(pos_deg, neg_deg, outdir)

    # --- 2. Degree ratios ---
    all_domains = set(pos_deg.keys()) | set(neg_deg.keys())
    degree_ratios = []
    for d in all_domains:
        n_pos = pos_deg.get(d, 0)
        n_neg = neg_deg.get(d, 0)
        total = n_pos + n_neg
        if total > 0:
            degree_ratios.append(n_pos / total)

    ratios_arr = np.array(degree_ratios)
    n_pure_positive = int(np.sum(ratios_arr == 1.0))
    n_pure_negative = int(np.sum(ratios_arr == 0.0))

    report["degree_ratios"] = {
        "n_total_domains": len(all_domains),
        "n_pure_positive": n_pure_positive,
        "frac_pure_positive": n_pure_positive / max(1, len(all_domains)),
        "n_pure_negative": n_pure_negative,
        "frac_pure_negative": n_pure_negative / max(1, len(all_domains)),
        "mean_ratio": float(np.mean(ratios_arr)) if len(ratios_arr) > 0 else 0.0,
        "median_ratio": float(np.median(ratios_arr)) if len(ratios_arr) > 0 else 0.0,
        "warning": (
            "High fraction of pure-positive/pure-negative domains enables "
            "degree-based shortcut predictions (Bernett et al. 2024)"
            if (n_pure_positive + n_pure_negative) / max(1, len(all_domains)) > 0.5
            else "OK"
        ),
    }

    plot_degree_ratios(degree_ratios, outdir)

    # --- 3. Hub / study bias ---
    total_deg = Counter()
    for d in all_domains:
        total_deg[d] = pos_deg.get(d, 0) + neg_deg.get(d, 0)

    # Get pfam_id for readable output
    pfam_lookup = dict(conn.execute("SELECT id, pfam_id FROM domain").fetchall())

    top_20 = total_deg.most_common(20)
    top_20_info = []
    for domain_id, total in top_20:
        top_20_info.append({
            "domain_id": domain_id,
            "pfam_id": pfam_lookup.get(domain_id, "?"),
            "degree_positive": pos_deg.get(domain_id, 0),
            "degree_negative": neg_deg.get(domain_id, 0),
            "degree_total": total,
        })

    gini = gini_coefficient(list(total_deg.values()))

    report["hub_study_bias"] = {
        "gini_coefficient": gini,
        "gini_interpretation": (
            "High inequality (>0.6) suggests hub-dominated network; "
            "models may learn hub identity rather than interaction features"
        ),
        "top_20_domains": top_20_info,
    }

    # --- Summary ---
    n_pos_ddis = conn.execute("SELECT COUNT(*) FROM domain_domain_interaction WHERE negative = 0").fetchone()[0]
    n_neg_ddis = conn.execute("SELECT COUNT(*) FROM domain_domain_interaction WHERE negative = 1").fetchone()[0]

    report["summary"] = {
        "n_positive_ddis": n_pos_ddis,
        "n_negative_ddis": n_neg_ddis,
        "pos_neg_ratio": n_pos_ddis / max(1, n_neg_ddis),
        "reference": "Bernett et al. (2024) Cracking the black box of deep PPI prediction. Brief Bioinform 25(2):bbae076",
    }

    conn.close()

    # Write JSON report
    report_path = outdir / "bias_report.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)

    # Write human-readable summary
    txt_path = outdir / "bias_report.txt"
    with open(txt_path, "w") as f:
        f.write("=== DDI Bias Analysis Report ===\n\n")
        f.write(f"Positive DDIs: {n_pos_ddis}\n")
        f.write(f"Negative DDIs: {n_neg_ddis}\n")
        f.write(f"Ratio: {n_pos_ddis / max(1, n_neg_ddis):.2f}\n\n")

        f.write("--- Degree Distribution ---\n")
        dd = report["degree_distribution"]
        f.write(f"Positive: mean={dd['pos_degree_mean']:.1f}, median={dd['pos_degree_median']:.1f} ({dd['n_domains_in_positive']} domains)\n")
        f.write(f"Negative: mean={dd['neg_degree_mean']:.1f}, median={dd['neg_degree_median']:.1f} ({dd['n_domains_in_negative']} domains)\n")
        f.write(f"KS test: statistic={dd['ks_statistic']:.4f}, p-value={dd['ks_p_value']:.4e}\n\n")

        f.write("--- Degree Ratios ---\n")
        dr = report["degree_ratios"]
        f.write(f"Pure positive domains: {dr['n_pure_positive']}/{dr['n_total_domains']} ({dr['frac_pure_positive']:.1%})\n")
        f.write(f"Pure negative domains: {dr['n_pure_negative']}/{dr['n_total_domains']} ({dr['frac_pure_negative']:.1%})\n")
        f.write(f"Status: {dr['warning']}\n\n")

        f.write("--- Hub / Study Bias ---\n")
        f.write(f"Gini coefficient: {report['hub_study_bias']['gini_coefficient']:.3f}\n")
        f.write(f"Top-5 domains by degree:\n")
        for d in report["hub_study_bias"]["top_20_domains"][:5]:
            f.write(f"  {d['pfam_id']}: total={d['degree_total']} (pos={d['degree_positive']}, neg={d['degree_negative']})\n")

    print(f"Bias report written to {outdir}")
    print(f"  {report_path}")
    print(f"  {txt_path}")
    print(f"  {outdir / 'degree_distribution.png'}")
    print(f"  {outdir / 'degree_ratios.png'}")


if __name__ == "__main__":
    main()
