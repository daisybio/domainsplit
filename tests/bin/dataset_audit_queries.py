#!/usr/bin/env python3
"""Read-only verification of the DATASET_AUDIT.md findings against a finished run.

Runs on the cluster, next to the published results. Opens every database with
``mode=ro`` and writes nothing anywhere except ``--out``; it cannot modify a run.

    python3 dataset_audit_queries.py \
        --results /nfs/proj/.../running_everything_now_fixed_with_tiers \
        --log     /nfs/proj/.../running_everything_now_fixed_with_tiers/.nextflow.log.everything.tiers \
        --out     audit_out

Then copy back ``audit_out/audit_report.txt`` (and the small TSVs beside it).

What each section decides, keyed to DATASET_AUDIT.md:

  F3  exact count of split rows whose CSV label contradicts the published
      ``domain_domain_interaction.negative`` -- the mislabel this audit predicts
  F1  is ``test_balanced``'s negative set a subset of ``test_realistic``'s, and
      are the positives identical
  F2  are the ``_hcni`` and non-``_hcni`` test sets identical, per ratio
  F4  external_test label composition by source; contamination count (expect 0);
      abundance of negatome vs single_domain_ppi families
  F5  true prevalence of the test domain universe
  F6  families shared between splits of one method (expect 0 for minimal_leakage*)
  F7  protein node/edge overlap between a method's splits
  F8  reversed rows vs unordered pairs in protein_protein_interaction, and the
      STRING score distribution
"""

import argparse
import csv
import os
import re
import sqlite3
import sys
from collections import Counter, defaultdict

MASTER = "domainsplit.sqlite3"
INTERNAL = ["random", "minimal_leakage", "minimal_leakage_hcni"]
RATIOS = ["test_balanced", "test_realistic"]


def ro(path):
    """Read-only connection. Fails rather than creating a database."""
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return sqlite3.connect(f"file:{os.path.abspath(path)}?mode=ro", uri=True)


def attach_ro(conn, path, name):
    conn.execute("ATTACH DATABASE ? AS " + name, (f"file:{os.path.abspath(path)}?mode=ro",))


class Report:
    def __init__(self, fh):
        self.fh = fh

    def head(self, text):
        self.fh.write(f"\n{'=' * 78}\n{text}\n{'=' * 78}\n")

    def line(self, text=""):
        self.fh.write(f"{text}\n")

    def rows(self, header, rows, limit=40):
        self.line("  " + " | ".join(header))
        for r in list(rows)[:limit]:
            self.line("  " + " | ".join("" if v is None else str(v) for v in r))


# --------------------------------------------------------------------------- F3

def canonical(a, b):
    """Pfam pair ordered by accession number -- mirrors bin/ddi_db_utils.py."""
    def key(p):
        digits = "".join(c for c in p if c.isdigit())
        return (0, int(digits)) if digits else (1, p)
    return tuple(sorted((a, b), key=key))


def ddi_labels(master):
    """{(pfam_a, pfam_b): (negative, source)} from the master."""
    out = {}
    for a, b, neg, src in master.execute(
        "SELECT da.pfam_id, db.pfam_id, ddi.negative, ddi.source "
        "FROM domain_domain_interaction AS ddi "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b"
    ):
        out[canonical(a, b)] = (int(neg), src)
    return out


def find_split_csvs(results):
    """[(row_id, filename, path)] for every published instance-level split CSV."""
    root = os.path.join(results, "ppi_splitting")
    found = []
    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in filenames:
            if fn.endswith("_instances.csv"):
                parts = os.path.relpath(dirpath, root).split(os.sep)
                found.append((parts[0] if parts else "?", fn, os.path.join(dirpath, fn)))
    return sorted(found)


def check_f3(rep, results, master, outdir):
    rep.head("F3 -- CSV label vs published domain_domain_interaction.negative")
    labels = ddi_labels(master)
    csvs = find_split_csvs(results)
    if not csvs:
        rep.line("  no *_instances.csv found under results/ppi_splitting -- skipping")
        return
    mismatch_rows = []
    per_file = Counter()
    for row_id, fn, path in csvs:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh)
            if not reader.fieldnames or "family1" not in reader.fieldnames:
                continue
            for r in reader:
                pair = canonical(r["family1"].strip(), r["family2"].strip())
                known = labels.get(pair)
                if known is None:
                    continue
                negative, source = known
                label = int(r["label"])
                if negative != (1 - label):
                    per_file[(row_id, fn)] += 1
                    if len(mismatch_rows) < 5000:
                        mismatch_rows.append((row_id, fn, pair[0], pair[1], label, negative, source))
    rep.line(f"  files scanned: {len(csvs)}")
    rep.line(f"  MISMATCHED ROWS TOTAL: {sum(per_file.values())}")
    if per_file:
        rep.line("  (a CSV negative published as a positive, or vice versa)")
        rep.rows(["row", "file", "n"], [(a, b, n) for (a, b), n in per_file.most_common()])
        out = os.path.join(outdir, "f3_label_mismatches.tsv")
        with open(out, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["row_id", "file", "pfam_a", "pfam_b", "csv_label", "db_negative", "db_source"])
            w.writerows(mismatch_rows)
        rep.line(f"  detail -> {out}")
    else:
        rep.line("  NO PROBLEM FOUND: every CSV label agrees with the published label.")


def check_f3_log(rep, log_path):
    rep.head("F3 -- INGEST_SAMPLED_NEGATIVES collision counters (from the run log)")
    if not log_path or not os.path.exists(log_path):
        rep.line("  no --log given or file missing; grep the work dirs instead:")
        rep.line("    grep -rh 'already_held_by_protected' $NXF_WORK/*/*/.command.log")
        return
    pat = re.compile(r"(already_held_by_protected|n_ddis_sampled_negative|no_ddi_row|unknown_instance|"
                     r"dropped_protected|conflict_pairs|sources_dropped)\s*=\s*(\d+)")
    hits = Counter()
    with open(log_path, errors="replace") as fh:
        for line in fh:
            m = pat.search(line)
            if m:
                hits[m.group(1)] = max(hits[m.group(1)], int(m.group(2)))
    if hits:
        rep.rows(["counter", "value"], sorted(hits.items()))
        if hits.get("already_held_by_protected"):
            rep.line("  ^ non-zero already_held_by_protected CONFIRMS the F3 mechanism on this run.")
    else:
        rep.line("  counters not present in this log (they are printed by the task, not the driver);")
        rep.line("  grep the work dirs: grep -rh 'already_held_by_protected' $NXF_WORK/*/*/.command.log")


# ----------------------------------------------------------------------- F1/F2

def pair_sets(db_path):
    """(positives, negatives) as {(pfam_a, pfam_b, instance_a, instance_b)}."""
    con = ro(db_path)
    q = ("SELECT da.pfam_id, dz.pfam_id, m.instance_id_a, m.instance_id_b, ddi.negative "
         "FROM ddi_split_membership AS m "
         "JOIN domain_domain_interaction AS ddi ON ddi.id = m.ddi_id "
         "JOIN domain AS da ON da.id = ddi.domain_id_a "
         "JOIN domain AS dz ON dz.id = ddi.domain_id_b")
    pos, neg = set(), set()
    for a, b, ia, ib, negative in con.execute(q):
        (neg if int(negative) else pos).add((a, b, ia, ib))
    con.close()
    return pos, neg


def check_f1(rep, results):
    rep.head("F1 -- is test_balanced a subset of test_realistic, and are positives equal")
    rep.line("  method | bal_neg | real_neg | bal_neg_not_in_real | bal_pos | real_pos | pos_equal")
    for method in INTERNAL:
        paths = [os.path.join(results, "databases", method, f"{s}.sqlite3") for s in RATIOS]
        if not all(os.path.exists(p) for p in paths):
            rep.line(f"  {method}: missing split DB(s) -- skipped")
            continue
        (bp, bn), (rp, rn) = pair_sets(paths[0]), pair_sets(paths[1])
        rep.line(f"  {method} | {len(bn)} | {len(rn)} | {len(bn - rn)} | "
                 f"{len(bp)} | {len(rp)} | {bp == rp}")
    rep.line("  bal_neg_not_in_real == 0 and pos_equal == True is what F1 requires.")


def check_f2(rep, results):
    rep.head("F2 -- are the _hcni and non-_hcni test sets identical")
    for split in RATIOS:
        a = os.path.join(results, "databases", "minimal_leakage", f"{split}.sqlite3")
        b = os.path.join(results, "databases", "minimal_leakage_hcni", f"{split}.sqlite3")
        if not (os.path.exists(a) and os.path.exists(b)):
            rep.line(f"  {split}: missing DB -- skipped")
            continue
        (ap, an), (bp, bn) = pair_sets(a), pair_sets(b)
        rep.line(f"  {split}: positives identical = {ap == bp}; negatives identical = {an == bn}; "
                 f"negatives differing = {len(an ^ bn)}")
    rep.line("  Expected today: test_realistic identical, test_balanced differing (the F2 bug).")


# --------------------------------------------------------------------------- F4

def check_f4(rep, master):
    rep.head("F4 -- external_test composition, contamination, abundance")

    rep.line("\n  (a) train/validation negatives that are external positives -- expect 0")
    n = master.execute(
        "SELECT COUNT(*) FROM ddi_split_membership m "
        "JOIN domain_domain_interaction d ON d.id = m.ddi_id "
        "WHERE m.split IN ('train','validation') AND d.negative = 1 "
        "AND d.source NOT LIKE '%sampled_negative%' AND d.source NOT LIKE '%3did%'"
    ).fetchone()[0]
    rep.line(f"      {n}")
    if n == 0:
        rep.line("      NO PROBLEM FOUND: the protected-source rule held.")

    rep.line("\n  (f/g) external_test composition by source and label")
    rep.rows(["negative", "source", "ddis", "instance_pairs"], master.execute(
        "SELECT d.negative, d.source, COUNT(DISTINCT d.id), COUNT(*) "
        "FROM ddi_split_membership m "
        "JOIN domain_domain_interaction d ON d.id = m.ddi_id "
        "WHERE m.method = 'external_test' AND m.split = 'test' "
        "GROUP BY 1, 2 ORDER BY 4 DESC"))

    rep.line("\n  (d) abundance: families of negatome vs single_domain_ppi")
    rep.rows(["class", "families", "mean_instances", "mean_3did_degree"], master.execute("""
        WITH src AS (
          SELECT DISTINCT da.pfam_id AS pfam,
                 CASE WHEN d.source LIKE '%negatome%' THEN 'negatome'
                      WHEN d.source LIKE '%single_domain_ppi%' THEN 'single_domain_ppi' END AS cls
          FROM domain_domain_interaction d
          JOIN domain da ON da.id IN (d.domain_id_a, d.domain_id_b)
          WHERE d.source LIKE '%negatome%' OR d.source LIKE '%single_domain_ppi%')
        SELECT cls, COUNT(*),
               AVG((SELECT COUNT(*) FROM domain_protein_map dm
                    JOIN domain dd ON dd.id = dm.domain_id WHERE dd.pfam_id = src.pfam)),
               AVG((SELECT COUNT(*) FROM domain_domain_interaction d2
                    JOIN domain d2a ON d2a.id IN (d2.domain_id_a, d2.domain_id_b)
                    WHERE d2a.pfam_id = src.pfam AND d2.source LIKE '%3did%'))
        FROM src WHERE cls IS NOT NULL GROUP BY 1"""))

    rep.line("\n  (c) domain coverage of the parent protein (single-domain assumption)")
    rep.rows(["coverage_bucket", "instances"], master.execute("""
        SELECT CAST(10.0 * (dm.end_pos - dm.start_pos + 1) / LENGTH(p.sequence) AS INT) / 10.0,
               COUNT(*)
        FROM domain_protein_map dm JOIN protein p ON p.id = dm.protein_id
        WHERE p.sequence IS NOT NULL AND LENGTH(p.sequence) > 0
        GROUP BY 1 ORDER BY 1"""))


# ------------------------------------------------------------------------ F5/F6

def check_f5(rep, results):
    rep.head("F5 -- true prevalence of the test domain universe")
    for method in INTERNAL:
        path = os.path.join(results, "databases", method, "test_realistic.sqlite3")
        if not os.path.exists(path):
            continue
        con = ro(path)
        domains = con.execute("SELECT COUNT(*) FROM domain").fetchone()[0]
        pos = con.execute(
            "SELECT COUNT(DISTINCT ddi_id) FROM ddi_split_membership m "
            "JOIN domain_domain_interaction d ON d.id = m.ddi_id WHERE d.negative = 0").fetchone()[0]
        con.close()
        possible = domains * (domains + 1) // 2
        rep.line(f"  {method}: {domains} domains -> {possible} possible pairs, {pos} positive DDIs "
                 f"= 1:{possible // pos if pos else 0}")


def check_f6(rep, master):
    rep.head("F6 -- families shared between splits of one method (expect 0 for minimal_leakage*)")
    rep.rows(["method", "split_a", "split_b", "negative", "shared_families"], master.execute("""
        WITH fam AS (
          SELECT DISTINCT m.method, m.split, da.pfam_id AS pfam, d.negative
          FROM ddi_split_membership m
          JOIN domain_domain_interaction d ON d.id = m.ddi_id
          JOIN domain da ON da.id IN (d.domain_id_a, d.domain_id_b))
        SELECT a.method, a.split, b.split, a.negative, COUNT(DISTINCT a.pfam)
        FROM fam a JOIN fam b
          ON a.method = b.method AND a.pfam = b.pfam AND a.split < b.split AND a.negative = b.negative
        GROUP BY 1,2,3,4 ORDER BY 1,2,3,4"""), limit=80)
    rep.line("  Rows for random are expected (edge split). Rows for minimal_leakage* are leakage.")

    rep.line("\n  clan cluster sizes (the actual disjointness unit)")
    rep.rows(["clan", "families"], master.execute(
        "SELECT COALESCE(clan, '(none)'), COUNT(DISTINCT domain_id) FROM domain_protein_map "
        "GROUP BY 1 ORDER BY 2 DESC"), limit=15)


# ------------------------------------------------------------------------ F7/F8

def check_f7(rep, results):
    rep.head("F7 -- protein node / PPI edge overlap between splits")
    for method in INTERNAL:
        d = os.path.join(results, "databases", method)
        tr, te = os.path.join(d, "train.sqlite3"), os.path.join(d, "test_balanced.sqlite3")
        if not (os.path.exists(tr) and os.path.exists(te)):
            continue
        con = ro(tr)
        attach_ro(con, te, "te")
        n_tr = con.execute("SELECT COUNT(*) FROM main.protein").fetchone()[0]
        n_te = con.execute("SELECT COUNT(*) FROM te.protein").fetchone()[0]
        shared = con.execute(
            "SELECT COUNT(*) FROM main.protein a JOIN te.protein b ON a.uniprot_id = b.uniprot_id"
        ).fetchone()[0]
        e_tr = con.execute("SELECT COUNT(*) FROM main.protein_protein_interaction").fetchone()[0]
        e_te = con.execute("SELECT COUNT(*) FROM te.protein_protein_interaction").fetchone()[0]
        con.close()
        rep.line(f"  {method}: train_nodes={n_tr} test_nodes={n_te} shared_nodes={shared} "
                 f"({100.0 * shared / max(1, min(n_tr, n_te)):.1f}% of the smaller) "
                 f"train_edges={e_tr} test_edges={e_te}")


def check_f8(rep, master):
    rep.head("F8 -- PPI duplication and score distribution")
    rows, reversed_rows = master.execute(
        "SELECT COUNT(*), SUM(CASE WHEN protein_id_a > protein_id_b THEN 1 ELSE 0 END) "
        "FROM protein_protein_interaction").fetchone()
    unordered = master.execute(
        "SELECT COUNT(*) FROM (SELECT DISTINCT MIN(protein_id_a, protein_id_b), "
        "MAX(protein_id_a, protein_id_b) FROM protein_protein_interaction)").fetchone()[0]
    proteins = master.execute("SELECT COUNT(*) FROM protein").fetchone()[0]
    rep.line(f"  rows={rows} reversed_rows={reversed_rows or 0} unordered_pairs={unordered} "
             f"proteins={proteins}")
    if unordered and rows / unordered > 1.5:
        rep.line(f"  ^ rows/unordered = {rows / unordered:.2f} -- CONFIRMS F8(ii): every edge stored twice.")
    else:
        rep.line("  NO PROBLEM FOUND for F8(ii): rows are not doubled.")
    rep.rows(["min", "avg", "max", "frac_below_400", "frac_below_700"], master.execute(
        "SELECT MIN(score), ROUND(AVG(score), 1), MAX(score), "
        "ROUND(SUM(score < 400) * 1.0 / COUNT(*), 3), ROUND(SUM(score < 700) * 1.0 / COUNT(*), 3) "
        "FROM protein_protein_interaction"))
    rep.line("  A large frac_below_400 means text-mined / predicted edges dominate (F8(i)).")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True, help="published results directory")
    ap.add_argument("--log", default=None, help="the run's .nextflow.log")
    ap.add_argument("--out", default="audit_out")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    master_path = os.path.join(args.results, MASTER)
    if not os.path.exists(master_path):
        sys.exit(f"{master_path} not found -- is --results the published output directory?")

    with open(os.path.join(args.out, "audit_report.txt"), "w") as fh:
        rep = Report(fh)
        rep.line(f"dataset audit verification -- results = {args.results}")
        master = ro(master_path)
        for fn, a in ((check_f3_log, (args.log,)), (check_f3, (args.results, master, args.out)),
                      (check_f1, (args.results,)), (check_f2, (args.results,)),
                      (check_f4, (master,)), (check_f5, (args.results,)), (check_f6, (master,)),
                      (check_f7, (args.results,)), (check_f8, (master,))):
            try:
                fn(rep, *a)
            except Exception as exc:  # one broken section must not lose the rest
                rep.line(f"  !! {fn.__name__} failed: {type(exc).__name__}: {exc}")
        master.close()
    print(f"wrote {os.path.join(args.out, 'audit_report.txt')}")


if __name__ == "__main__":
    main()
