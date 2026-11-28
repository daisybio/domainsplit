#!/usr/bin/env python3
"""AF3 prediction queue: every protein pair still missing a model, in level order.

Reads the master right after ENRICH_STRUCTURAL_AF (raw has_structure, before the
coverage-balance guard) and the af_metadata.csv it was sliced from.

Levels are fixed (LEVELS). Every (method, split) of the selected methods falls in
exactly one level; one that matches none fails, so a new split cannot drop out of
the queue silently. A protein pair is listed once, at the first level that needs
it, so level N lists only what it adds to levels 1..N-1. Unit: unordered protein
pair (one AF3 job), uniprot_a <= uniprot_b (equal = a homodimer job for two domains
of one protein). "Missing a model" = the pair has no row
in af_metadata.csv; pairs that have a model nobody could slice go to
unusable_models.tsv instead (re-predicting them is a separate decision).

Within a level, rows are ordered by CRC32(seed:uniprot_a:uniprot_b): a partly
submitted level is label-mixed, and a pair keeps its place across dataset versions.

Outputs (in --out_dir):
  queue.tsv            plain TSV, one header line, rows in level order.
                       Levels 1-2 only: awk -F'\t' 'NR == 1 || $1 <= 2'
  queue_summary.txt    one line per level with the pairs it adds
  unusable_models.tsv  pairs with a model but no sliced instance pair, with the
                       gap report's reason and the model's source/path
and --counts (af3_queue.json): the per-level counts, embedded into run_manifest.json.
"""
import argparse
import csv
import json
import os
import sqlite3
import sys
import zlib
from collections import defaultdict

# (level_name, method or None = every selected method, splits)
LEVELS = [
    ("minimal_leakage", "minimal_leakage", ("train", "validation", "test_matched_1_1")),
    ("minimal_leakage_hrni", "minimal_leakage_hrni", ("train", "validation", "test_matched_1_1")),
    ("random_ilp", "random_ilp", ("train", "validation", "test_matched_1_1")),
    ("random_hrni", "random_hrni", ("train", "validation", "test_matched_1_1")),
    ("external_test", "external_test", ("train", "validation", "test")),
    ("external_test_hrni", "external_test_hrni", ("train", "validation", "test")),
    ("random", "random", ("train", "validation", "test_uniform_1_1")),
    ("test_uniform_1_10", None, ("test_uniform_1_10",)),
]
COLUMNS = ["level", "level_name", "uniprot_a", "uniprot_b", "labels", "splits", "n_instance_pairs"]


def level_of(method, split):
    for i, (_name, m, splits) in enumerate(LEVELS, start=1):
        if (m is None or m == method) and split in splits:
            return i
    return None


def protein(instance_id):
    # <pfam_id>_<uniprot_id>_<start>_<end>
    return instance_id.split("_")[1]


def rank(seed, a, b):
    return zlib.crc32(f"{seed}:{a}:{b}".encode()), a, b


def load_af_pairs(path):
    """{(a, b) sorted: (source, path)} from af_metadata.csv."""
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            a, b = sorted((row["uniprot_id_a"], row["uniprot_id_b"]))
            out[(a, b)] = (row.get("source", ""), row.get("path", ""))
    return out


def load_reasons(path):
    """{(a, b): reason} from af3_to_model/all_methods.tsv (already sorted pairs)."""
    with open(path, newline="") as fh:
        return {(r["uniprot_id_a"], r["uniprot_id_b"]): r["reason"]
                for r in csv.DictReader(fh, delimiter="\t")}


def collect(db, methods):
    """Per protein pair: first level, and at that level its labels, splits and
    instance pairs; plus whether any instance pair (any level) has a structure."""
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    placeholders = ",".join("?" for _ in methods)
    rows = conn.execute(
        "SELECT m.method, m.split, m.instance_id_a, m.instance_id_b, m.has_structure, d.negative "
        "FROM ddi_split_membership m JOIN domain_domain_interaction d ON d.id = m.ddi_id "
        "WHERE m.instance_id_a IS NOT NULL AND m.instance_id_b IS NOT NULL "
        f"AND m.method IN ({placeholders})", methods).fetchall()
    conn.close()

    unmatched = sorted({f"{m}/{s}" for m, s, *_ in rows if level_of(m, s) is None})
    if unmatched:
        sys.exit(f"build_af3_queue: split(s) {', '.join(unmatched)} match no queue level; "
                 "add them to LEVELS in bin/build_af3_queue.py")

    first = {}                                  # pair -> first level
    structured = set()                          # pairs with >= 1 sliced instance pair
    for method, split, ia, ib, has, _neg in rows:
        pair = tuple(sorted((protein(ia), protein(ib))))
        lvl = level_of(method, split)
        first[pair] = min(first.get(pair, lvl), lvl)
        if has:
            structured.add(pair)
    detail = defaultdict(lambda: (set(), set(), set()))   # pair -> labels, splits, instance pairs
    for method, split, ia, ib, _has, neg in rows:
        pair = tuple(sorted((protein(ia), protein(ib))))
        if level_of(method, split) == first[pair]:
            labels, splits, inst = detail[pair]
            labels.add("neg" if neg else "pos")
            splits.add(f"{method}/{split}")
            inst.add((ia, ib))
    return first, detail, structured


def fmt_counts(c):
    s = f"{c['pos']:,} pos, {c['neg']:,} neg"
    return s + (f", {c['pairs']:,} pairs" if c["pairs"] != c["pos"] + c["neg"] else "")


def row_for(lvl, pair, detail):
    labels, splits, inst = detail[pair]
    return [lvl, LEVELS[lvl - 1][0], pair[0], pair[1],
            ",".join(sorted(labels, reverse=True)),   # "pos,neg"
            ",".join(sorted(splits)), len(inst)]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--db", required=True)
    p.add_argument("--af_metadata", required=True)
    p.add_argument("--gap_report", required=True, help="af3_to_model/all_methods.tsv")
    p.add_argument("--methods", required=True, help="comma-separated structure methods")
    p.add_argument("--seed", required=True)
    p.add_argument("--out_dir", required=True)
    p.add_argument("--counts", required=True, help="af3_queue.json (per-level counts)")
    args = p.parse_args(argv)
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]

    af_pairs = load_af_pairs(args.af_metadata)
    reasons = load_reasons(args.gap_report)
    first, detail, structured = collect(args.db, methods)

    by_level = defaultdict(list)
    unusable = []
    for pair, lvl in first.items():
        if pair not in af_pairs:
            by_level[lvl].append(pair)
        elif pair not in structured:
            unusable.append((lvl, pair))

    os.makedirs(args.out_dir, exist_ok=True)
    counts = []
    with open(os.path.join(args.out_dir, "queue.tsv"), "w") as fh:
        fh.write("\t".join(COLUMNS) + "\n")
        for lvl, (name, method, _splits) in enumerate(LEVELS, start=1):
            pairs = sorted(by_level[lvl], key=lambda pr: rank(args.seed, *pr))
            c = {"level": lvl, "name": name,
                 "selected": method is None or method in methods,
                 "pos": sum("pos" in detail[pr][0] for pr in pairs),
                 "neg": sum("neg" in detail[pr][0] for pr in pairs),
                 "pairs": len(pairs)}
            counts.append(c)
            for pr in pairs:
                fh.write("\t".join(map(str, row_for(lvl, pr, detail))) + "\n")

    with open(os.path.join(args.out_dir, "queue_summary.txt"), "w") as fh:
        selected = [c for c in counts if c["selected"]]
        for i, c in enumerate(selected):
            fh.write(f"{'+ ' if i else ''}{c['name']} ({fmt_counts(c)})\n")

    with open(os.path.join(args.out_dir, "unusable_models.tsv"), "w") as fh:
        fh.write("\t".join(COLUMNS + ["reason", "source", "path"]) + "\n")
        for lvl, pair in sorted(unusable):
            if pair not in reasons:
                sys.exit(f"build_af3_queue: {pair} has a model and no sliced instance pair, "
                         "but no row in the gap report")
            fh.write("\t".join(map(str, row_for(lvl, pair, detail) +
                                   [reasons[pair], *af_pairs[pair]])) + "\n")

    with open(args.counts, "w") as fh:
        json.dump({"seed": args.seed, "levels": counts, "unusable_models": len(unusable)},
                  fh, indent=2)
        fh.write("\n")
    total = sum(c["pairs"] for c in counts)
    print(f"build_af3_queue: {total} protein pairs queued over {len(LEVELS)} levels, "
          f"{len(unusable)} with an unusable model", flush=True)


if __name__ == "__main__":
    main()
