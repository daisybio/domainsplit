#!/usr/bin/env python3
"""Split QC on the final master: overlap between splits, test-variant overlap,
the external test set's hub profile and (optionally) esm3_structure coverage.

Unit of everything here is ``ddi_split_membership`` (method, split, ddi_id,
instance_id_a, instance_id_b); family / protein come from the instance id
(``{pfam}_{uniprot}_{start}_{end}``), the label from
``domain_domain_interaction.negative`` (same rule as write_run_manifest.py).

Outputs (all TSV, deterministically sorted):

* ``split_overlap.tsv``: per method and unordered split pair, the size of each
  split's set and the shared count, for families, positive / negative family
  pairs, proteins and instances.
* ``family_leaks.tsv``: for the ``--ilp-methods``, every family found in more
  than one of train / validation / the test_* splits (one partition). The ILP split guarantees there are
  none, so any row makes the script exit 1 (after writing every other file).
  The external ``test`` split is not part of that partition and is not checked.
* ``test_variant_overlap.tsv``: per method with internal tests, the 1:1 test
  vs ``test_uniform_1_10``, per label, at DDI-row and family-pair granularity.
* ``external_test_profile.tsv``: per ``--external-methods`` method, its ``test``
  split per (label, source token): sizes, family overlap with the other label
  and with train, and the train-positive degree of the pair's families.
* with ``--esm3-h5`` + ``--afdb-report``: ``esm3_structure_coverage.tsv``,
  ``esm3_structure_missing_instances.tsv``,
  ``esm3_structure_missing_accessions.tsv``. Missing coverage only warns.
"""

import argparse
import csv
import itertools
import os
import sqlite3
import statistics
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ddi_db_utils import split_sources  # noqa: E402
from utils_struct import extract_data_from_instance_id  # noqa: E402

# BUILD_EXTERNAL_TEST's split name (workflows/domainsplit.nf passes it 'test').
EXTERNAL_SPLIT = "test"
ONE_TO_ONE = ("test_matched_1_1", "test_uniform_1_1")
ONE_TO_TEN = "test_uniform_1_10"
LABELS = ("positive", "negative")


def log(msg):
    print(f"[split_qc] {msg}", file=sys.stderr, flush=True)


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(sorted(rows))


def partition(split):
    """The ILP partition a split belongs to, or None when the ILP did not make it.

    Every test_* split is drawn from the one test partition (the 1:1 and 1:10
    sets share their positives), so they may share families with each other.
    """
    if split in ("train", "validation"):
        return split
    return "test" if split.startswith("test_") else None


def parse(instance_id, cache={}):
    if instance_id not in cache:
        pfam, acc, _s, _e = extract_data_from_instance_id(instance_id)
        cache[instance_id] = (pfam, acc)
    return cache[instance_id]


def fam_pair(a, b):
    return (a, b) if a <= b else (b, a)


def membership(db):
    """Yield ``(method, rows)`` with rows ``(split, ddi_id, inst_a, inst_b, label, source)``."""
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    cur = conn.execute(
        "SELECT m.method, m.split, m.ddi_id, m.instance_id_a, m.instance_id_b, d.negative, d.source "
        "FROM ddi_split_membership m JOIN domain_domain_interaction d ON d.id = m.ddi_id "
        "ORDER BY m.method"
    )
    for method, group in itertools.groupby(cur, key=lambda r: r[0]):
        yield method, [(s, ddi, a, b, "negative" if neg else "positive", src)
                       for _m, s, ddi, a, b, neg, src in group]
    conn.close()


def split_sets(rows):
    """Per split: families, positive/negative family pairs, proteins, instances."""
    sets = defaultdict(lambda: defaultdict(set))
    for split, _ddi, a, b, label, _src in rows:
        (pa, ua), (pb, ub) = parse(a), parse(b)
        s = sets[split]
        s["family"].update((pa, pb))
        s[f"family_pair_{label}"].add(fam_pair(pa, pb))
        s["protein"].update((ua, ub))
        s["instance"].update((a, b))
    return sets


def overlap_rows(method, sets):
    units = ("family", "family_pair_positive", "family_pair_negative", "protein", "instance")
    for sa, sb in itertools.combinations(sorted(sets), 2):
        for unit in units:
            A, B = sets[sa][unit], sets[sb][unit]
            yield (method, sa, sb, unit, len(A), len(B), len(A & B))


def leak_rows(method, sets):
    where = defaultdict(set)
    for split, s in sets.items():
        if partition(split):
            for fam in s["family"]:
                where[fam].add(split)
    return [(method, fam, ",".join(sorted(sp))) for fam, sp in where.items()
            if len({partition(x) for x in sp}) > 1]


def variant_rows(method, rows):
    present = {r[0] for r in rows}
    one = [s for s in ONE_TO_ONE if s in present]
    if ONE_TO_TEN not in present or not one:
        return []
    ddi_rows = defaultdict(set)
    pairs = defaultdict(set)
    for split, ddi, a, b, label, _src in rows:
        if split in (one[0], ONE_TO_TEN):
            ddi_rows[split, label].add((ddi,) + fam_pair(a, b))
            pairs[split, label].add(fam_pair(parse(a)[0], parse(b)[0]))
    out = []
    for label in LABELS:
        for gran, d in (("ddi_row", ddi_rows), ("family_pair", pairs)):
            A, B = d[one[0], label], d[ONE_TO_TEN, label]
            out.append((method, one[0], label, gran, len(A), len(B), len(A & B)))
    return out


def degree_stats(values):
    if not values:
        return ("", "", "")
    return (min(values), statistics.median(values), max(values))


def external_rows(method, rows):
    train_pos_pairs, train_pairs, train_fams = set(), set(), set()
    classes = defaultdict(lambda: {"pairs": set(), "n_rows": 0})
    label_fams = defaultdict(set)
    for split, _ddi, a, b, label, src in rows:
        fp = fam_pair(parse(a)[0], parse(b)[0])
        if split == "train":
            train_pairs.add(fp)
            train_fams.update(fp)
            if label == "positive":
                train_pos_pairs.add(fp)
        elif split == EXTERNAL_SPLIT:
            label_fams[label].update(fp)
            for token in split_sources(src) + ["__total__"]:
                c = classes[label, token]
                c["pairs"].add(fp)
                c["n_rows"] += 1
    if not classes:
        # `--external_test_balance none` can legitimately end with no external rows;
        # that must not cost the internal overlap and leak reports.
        log(f"WARNING: external method {method!r} has no '{EXTERNAL_SPLIT}' membership rows; "
            "its external_test_profile.tsv rows are empty")
    degree = Counter()
    for fa, fb in train_pos_pairs:
        degree[fa] += 1
        if fb != fa:
            degree[fb] += 1
    out = []
    for (label, token), c in classes.items():
        fams = {f for p in c["pairs"] for f in p}
        other = label_fams["negative" if label == "positive" else "positive"]
        out.append((
            method, label, token, len(c["pairs"]), c["n_rows"], len(fams),
            len(fams & other), len(fams & train_fams), len(c["pairs"] & train_pairs),
            *degree_stats([degree[f] for f in fams]),
            statistics.median(min(degree[a], degree[b]) for a, b in c["pairs"]),
        ))
    return out


def h5_instances(path):
    import h5py

    with h5py.File(path, "r") as h5:
        return {inst for fam in h5 for inst in h5[fam]}


def read_afdb_report(path):
    with open(path, newline="") as f:
        return {r["uniprot_id"]: r["status"] for r in csv.DictReader(f, delimiter="\t")}


class Esm3Coverage:
    def __init__(self, h5_path, report_path):
        self.present = h5_instances(h5_path)
        self.status = read_afdb_report(report_path)
        self.cov = Counter()
        self.missing_rows = Counter()

    def add(self, method, rows):
        for split, _ddi, a, b, label, _src in rows:
            ok_a, ok_b = a in self.present, b in self.present
            self.cov[method, split, label, "total"] += 1
            self.cov[method, split, label, "both"] += ok_a and ok_b
            for inst in {x for x, ok in ((a, ok_a), (b, ok_b)) if not ok}:
                self.missing_rows[inst] += 1

    def reason(self, acc):
        status = self.status.get(acc)
        if status is None:
            return "not_in_afdb_report"
        # The protein's model is fine, so run_embeddings.py left the domain out
        # for its structure coverage (--min-coverage) or its length (--max-len).
        return "below_min_coverage" if status == "ok" else status

    def write(self, outdir):
        keys = {k[:3] for k in self.cov}
        write_tsv(os.path.join(outdir, "esm3_structure_coverage.tsv"),
                  ["method", "split", "label", "instance_pairs", "both_in_h5", "fraction"],
                  [(*k, self.cov[(*k, "total")], self.cov[(*k, "both")],
                    round(self.cov[(*k, "both")] / self.cov[(*k, "total")], 4)) for k in keys])
        inst_rows, by_acc = [], Counter()
        for inst, n in self.missing_rows.items():
            acc = parse(inst)[1]
            reason = self.reason(acc)
            inst_rows.append((inst, acc, reason, n))
            by_acc[acc, reason] += 1
        write_tsv(os.path.join(outdir, "esm3_structure_missing_instances.tsv"),
                  ["instance_id", "uniprot_id", "reason", "n_membership_rows"], inst_rows)
        write_tsv(os.path.join(outdir, "esm3_structure_missing_accessions.tsv"),
                  ["accession", "reason", "n_instances"], [(*k, n) for k, n in by_acc.items()])
        if inst_rows:
            reasons = Counter(r[2] for r in inst_rows)
            log(f"WARNING: {len(inst_rows)} membership instances have no esm3_structure embedding "
                f"({len(by_acc)} accessions): {dict(sorted(reasons.items()))}")


def tokens(s):
    return [t for t in (x.strip() for x in s.split(",")) if t]


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db", required=True)
    p.add_argument("--ilp-methods", default="", help="comma-separated methods split by the ILP")
    p.add_argument("--external-methods", default="", help="comma-separated methods with an external test split")
    p.add_argument("--esm3-h5", default="")
    p.add_argument("--afdb-report", default="")
    p.add_argument("--outdir", default=".")
    p.add_argument("--versions")
    p.add_argument("--process-name")
    args = p.parse_args(argv)

    if bool(args.esm3_h5) != bool(args.afdb_report):
        raise SystemExit("--esm3-h5 and --afdb-report come together")
    ilp, external = set(tokens(args.ilp_methods)), set(tokens(args.external_methods))
    esm3 = Esm3Coverage(args.esm3_h5, args.afdb_report) if args.esm3_h5 else None

    overlap, leaks, variants, profile, seen = [], [], [], [], set()
    for method, rows in membership(args.db):
        seen.add(method)
        sets = split_sets(rows)
        overlap.extend(overlap_rows(method, sets))
        if method in ilp:
            leaks.extend(leak_rows(method, sets))
        variants.extend(variant_rows(method, rows))
        if method in external:
            profile.extend(external_rows(method, rows))
        if esm3:
            esm3.add(method, rows)

    absent = (ilp | external) - seen
    if absent:
        raise SystemExit(f"methods with no membership rows in {args.db}: {', '.join(sorted(absent))}")

    o = args.outdir
    write_tsv(os.path.join(o, "split_overlap.tsv"),
              ["method", "split_a", "split_b", "unit", "n_a", "n_b", "n_shared"], overlap)
    write_tsv(os.path.join(o, "family_leaks.tsv"), ["method", "pfam_id", "splits"], leaks)
    write_tsv(os.path.join(o, "test_variant_overlap.tsv"),
              ["method", "split_1_1", "label", "granularity", "n_1_1", "n_1_10", "n_shared"], variants)
    write_tsv(os.path.join(o, "external_test_profile.tsv"),
              ["method", "label", "source", "family_pairs", "instance_pairs", "families",
               "families_shared_with_other_label", "families_in_train", "family_pairs_in_train",
               "train_pos_degree_min", "train_pos_degree_median", "train_pos_degree_max",
               "pair_min_train_pos_degree_median"], profile)
    if esm3:
        esm3.write(o)

    if args.versions:
        with open(args.versions, "w") as f:
            f.write(f'"{args.process_name}":\n')
            f.write(f"    python: {sys.version.split()[0]}\n")
            f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")

    if leaks:
        # The task fails, so family_leaks.tsv is not published; a failed run keeps its
        # work directory, which is where the full list stays.
        log(f"ERROR: ILP family-disjointness violated: {len(leaks)} (method, family) leaks; "
            "full list in family_leaks.tsv in this task's work directory")
        for row in sorted(leaks)[:50]:
            log("  " + "\t".join(row))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
