#!/usr/bin/env python3
"""Coverage-balance guard: structure coverage must not carry the label (F11).

Per (method, split) of every method (not only structure_data_methods:
ENRICH_STRUCTURAL_AF flags an instance-pair triple in every method that shares
it, so an unselected method can carry structures too), a DDI is covered when at least one
of its membership rows has has_structure = 1. With the guard on:

  * any of the DDI coverage fractions (both labels together, positives,
    negatives) below --min_coverage: every structure of that split is masked
    (action `masked_all`); the split then carries no structure feature at all;
  * otherwise the label with the higher fraction is masked down to the lower
    one (action `balanced`): it keeps floor(f_low * n_high) of its covered DDIs,
    chosen by CRC32(seed:pfam_a:pfam_b) rank, so the choice does not depend on
    the surrogate ddi id and is the same in every split of a run.

A label without any DDI in the split is left out of the checks (nothing to
compare it with). floor() keeps the masked label at or below the other's
fraction; with few DDIs that can undershoot it (1/1 vs 3/10 keeps 0 of 1).
Masking a DDI sets has_structure = 0 and z_score = NULL on all
its rows in that (method, split) and deletes its ddi_interaction_confirmed row
for that method (one row per method, so a DDI masked in any split loses it).
has_structure_raw keeps the unmasked state on every row; structures.h5 is left
untouched. With the guard off nothing is masked, but the report is still
written, so the imbalance stays visible.

Runs after ANNOTATE_DDI: the contact potential still uses all train positives
(their labels are known legitimately), which does not decide which rows carry
a feature.
"""
import argparse
import json
import shutil
import sqlite3
import sys
import zlib
from collections import defaultdict

COLUMNS = ["method", "split", "label", "ddis", "covered_before", "covered_after", "masked",
           "fraction_before", "fraction_after", "action"]


def frac(num, den):
    return num / den if den else float("nan")


def plan_split(ddis, enabled, min_coverage):
    """ddis: {ddi_id: (label, covered, rank_key)} -> (action, set of ddi_ids to mask)."""
    n, c = defaultdict(int), defaultdict(int)
    for label, covered, _key in ddis.values():
        n[label] += 1
        c[label] += covered
    present = [lab for lab in ("positive", "negative") if n[lab]]
    if not enabled:
        return "off", set()
    total = frac(sum(c.values()), sum(n.values()))
    if total < min_coverage or any(frac(c[lab], n[lab]) < min_coverage for lab in present):
        return "masked_all", {d for d, (_l, cov, _k) in ddis.items() if cov}
    if len(present) < 2:
        return "single_label", set()
    # Integer comparison c_pos/n_pos vs c_neg/n_neg, no float ties.
    hi, lo = (("positive", "negative") if c["positive"] * n["negative"] > c["negative"] * n["positive"]
              else ("negative", "positive"))
    keep = c[lo] * n[hi] // n[lo]
    covered_hi = sorted((key, d) for d, (lab, cov, key) in ddis.items() if lab == hi and cov)
    return "balanced", {d for _key, d in covered_hi[keep:]}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--db_in", required=True)
    p.add_argument("--db_out", required=True)
    p.add_argument("--seed", required=True)
    p.add_argument("--enabled", required=True, choices=["true", "false"])
    p.add_argument("--min_coverage", type=float, required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--versions")
    p.add_argument("--process_name", default="STRUCTURE_COVERAGE_BALANCE")
    args = p.parse_args(argv)
    enabled = args.enabled == "true"
    if not 0 <= args.min_coverage <= 1:
        sys.exit(f"--min_coverage must be in [0, 1], got {args.min_coverage}")

    shutil.copy(args.db_in, args.db_out)
    conn = sqlite3.connect(args.db_out)
    cols = {r[1] for r in conn.execute("PRAGMA table_info(ddi_split_membership)")}
    if "has_structure_raw" not in cols:
        sys.exit("ddi_split_membership has no has_structure_raw column (INIT_DOMAINSPLIT_DB schema)")
    conn.execute("UPDATE ddi_split_membership SET has_structure_raw = has_structure")

    rows = conn.execute(
        "SELECT m.method, m.split, m.ddi_id, d.negative, MAX(m.has_structure), da.pfam_id, db.pfam_id "
        "FROM ddi_split_membership m "
        "JOIN domain_domain_interaction d ON d.id = m.ddi_id "
        "JOIN domain da ON da.id = d.domain_id_a JOIN domain db ON db.id = d.domain_id_b "
        "GROUP BY m.method, m.split, m.ddi_id").fetchall()
    splits = defaultdict(dict)
    for method, split, ddi_id, neg, covered, pfam_a, pfam_b in rows:
        key = zlib.crc32(f"{args.seed}:{pfam_a}:{pfam_b}".encode())
        splits[(method, split)][ddi_id] = ("negative" if neg else "positive", int(covered), (key, pfam_a, pfam_b))

    report, actions, masked_confirmed = [], {}, set()
    total_rows = 0
    for (method, split), ddis in sorted(splits.items()):
        action, mask = plan_split(ddis, enabled, args.min_coverage)
        actions[f"{method}/{split}"] = action
        for ddi_id in mask:
            total_rows += conn.execute(
                "UPDATE ddi_split_membership SET has_structure = 0, z_score = NULL "
                "WHERE method = ? AND split = ? AND ddi_id = ? AND has_structure = 1",
                (method, split, ddi_id)).rowcount
            masked_confirmed.add((ddi_id, method))
        for label in ("positive", "negative"):
            sel = [(d, cov) for d, (lab, cov, _k) in ddis.items() if lab == label]
            before = sum(cov for _d, cov in sel)
            masked = sum(1 for d, cov in sel if cov and d in mask)
            report.append([method, split, label, len(sel), before, before - masked, masked,
                           round(frac(before, len(sel)), 4), round(frac(before - masked, len(sel)), 4),
                           action])
    conn.executemany("DELETE FROM ddi_interaction_confirmed WHERE ddi_id = ? AND method = ?",
                     sorted(masked_confirmed))
    conn.commit()
    conn.close()

    with open(args.report, "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for r in report:
            fh.write("\t".join(map(str, r)) + "\n")
    summary = {
        "enabled": enabled,
        "min_coverage": args.min_coverage,
        "masked_ddis": sum(r[6] for r in report),
        "masked_rows": total_rows,
        "actions": actions,
    }
    with open(args.summary, "w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print(f"balance_structure_coverage: enabled={enabled}, {summary['masked_ddis']} DDI(s) masked "
          f"over {len(splits)} split(s), {total_rows} row(s)", flush=True)

    if args.versions:
        with open(args.versions, "w") as f:
            f.write(f'"{args.process_name}":\n')
            f.write(f"    python: {sys.version.split()[0]}\n")
            f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
