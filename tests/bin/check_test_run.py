#!/usr/bin/env python3
"""Assert the invariants of a finished `-profile test` run against its outputs.

The nf-test file calls this and fails on a non-zero exit, so the assertions live
here rather than in Groovy: they are about *database contents*, which is what a
`-stub` path snapshot could never check, and they are the properties the
integration is supposed to guarantee.

Checked:

1. the expected `databases/<method>/<split>.sqlite3` set exists, and nothing else;
2. every split DB is a subset of the master and holds only its own split's DDIs;
3. `3did` never shares a source list with an external source (the protected-source
   drop rule) and no DDI is both positive and negative;
4. at least one external DDI kept several source labels (`PPIDM,PPIDM_Gold`), so
   the merge path is exercised, or the run is reported as not covering it;
5. every DDI in the master DB has instances for both families (the prune invariant);
6. the external test split holds only external DDIs, each represented by instance
   pairs whose parent proteins differ;
7. the reports exist, and reports/run_manifest.json's per-split counts match
   the master's membership;
8. every published embedding HDF5 satisfies the benchmark's own key contract --
   `h5[pfam_id][COALESCE(instance_id, 'r' || rowid)]` resolves for the
   master's domain instances -- and carries its root attributes. This is the seam
   that yields zero training rows *silently* if it drifts, so it is checked here
   as well as in tests/python/test_export_domain_embeddings.py;
9. column types survived insertion: `protein_protein_interaction.score` is REAL
   and `domain_protein_map.start_pos/end_pos` are INTEGER in every database (a
   TEXT score once crashed domainbenchmark's confidence filter mid-run);
10. `minimal_leakage*` shares no Pfam family between any two of its
    train / validation / test splits (DATASET_AUDIT.md F6);
11. structure (only when the run sliced any): every `has_structure = 1` triple
    resolves in structures/structures.h5 and the file holds nothing else; a
    z-score only sits on a row with a structure and is identical across every
    split row of its (method, instance pair); `ddi_interaction_confirmed` only
    votes on DDIs its method scored; each published scoring matrix is a
    symmetric 20x20 with frequencies summing to 1 and T_DB consistent; split
    databases carry only their own method's confirmation rows.

Usage: ``tests/bin/check_test_run.py --outdir results``
"""

import argparse
import csv
import json
import os
import sqlite3
import sys

# `params.ppi_splitting_multi_negset = true` (the default) asks ppi-splitting for
# both ILP negative sets per ILP row, so each of those rows yields a plain
# directory and an `_hcni` one; `random` additionally gets `random_ilp` and
# `random_hcni` next to its uniform baseline: 7 method directories, 26 split
# databases. Set `--ppi_splitting_multi_negset false` and this drops back to 3 / 11.
# Published names, i.e. after publishedSplit() in workflows/domainsplit.nf: the
# 1:1 test set is named for its negative set's sampler -- uniform for `random`,
# the bias-matched ILP for the ILP rows -- and the 1:10 one is uniform
# everywhere (DATASET_AUDIT.md F1).
UNIFORM_SPLITS = ["train", "validation", "test_uniform_1_1", "test_uniform_1_10"]
MATCHED_SPLITS = ["train", "validation", "test_matched_1_1", "test_uniform_1_10"]
EXTERNAL_SPLITS = ["train", "validation", "test"]
EXPECTED = {
    "random": UNIFORM_SPLITS,
    "random_ilp": MATCHED_SPLITS,
    "random_hcni": MATCHED_SPLITS,
    "minimal_leakage": MATCHED_SPLITS,
    "minimal_leakage_hcni": MATCHED_SPLITS,
    "external_test": EXTERNAL_SPLITS,
    "external_test_hcni": EXTERNAL_SPLITS,
}
EXTERNAL_SOURCES = ("single_domain_ppi", "PPIDM", "negatome")
# Root attributes export_domain_embeddings.py writes. There is no run identifier:
# the files key on the Pfam accession, which needs no cross-run guard.
EMBEDDING_ATTRS = ("model", "pooling", "dim", "dtype", "key_layout",
                   "n_domains", "n_instances")

failures = []
notes = []


def check(condition, message):
    if not condition:
        failures.append(message)


def sources_of(row):
    return {s.strip() for s in (row or "").split(",") if s.strip()}


def ddi_rows(conn):
    return conn.execute(
        "SELECT da.pfam_id, db.pfam_id, ddi.negative, ddi.source "
        "FROM domain_domain_interaction AS ddi "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b"
    ).fetchall()


def check_layout(outdir):
    root = os.path.join(outdir, "databases")
    check(os.path.isdir(root), f"{root} does not exist")
    if not os.path.isdir(root):
        return {}
    found = {}
    for method in sorted(os.listdir(root)):
        splits = sorted(
            f[: -len(".sqlite3")] for f in os.listdir(os.path.join(root, method))
            if f.endswith(".sqlite3")
        )
        found[method] = splits
    check(set(found) == set(EXPECTED), f"method directories are {sorted(found)}, expected {sorted(EXPECTED)}")
    for method, splits in EXPECTED.items():
        check(found.get(method) == sorted(splits),
              f"databases/{method} holds {found.get(method)}, expected {sorted(splits)}")
    return found


def check_master(path):
    conn = sqlite3.connect(path)
    rows = ddi_rows(conn)
    check(bool(rows), "the master database has no DDIs at all")

    seen_labels = {}
    merged = 0
    for pfam_a, pfam_b, negative, source in rows:
        sources = sources_of(source)
        check(sources, f"{pfam_a}/{pfam_b} has an empty source list")
        protected = {s for s in sources if s in ("3did", "sampled_negative")}
        external = {s for s in sources if s.split("_")[0] in ("PPIDM",) or s in EXTERNAL_SOURCES}
        check(not (protected and external),
              f"{pfam_a}/{pfam_b} mixes protected {sorted(protected)} with external {sorted(external)}")
        if len(sources) > 1:
            merged += 1
        key = tuple(sorted((pfam_a, pfam_b)))
        check(key not in seen_labels or seen_labels[key] == negative,
              f"{pfam_a}/{pfam_b} appears with both labels")
        seen_labels[key] = negative

    if merged == 0:
        notes.append("no DDI kept more than one source label, so the merge path is untested "
                     "by this fixture (expected when PPIDM predicts each pair under one class)")

    # The prune invariant: every family of every surviving DDI has an instance.
    orphans = conn.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction AS ddi WHERE "
        "  ddi.domain_id_a NOT IN (SELECT domain_id FROM domain_protein_map) OR "
        "  ddi.domain_id_b NOT IN (SELECT domain_id FROM domain_protein_map)"
    ).fetchone()[0]
    check(orphans == 0, f"{orphans} DDIs survived with a family that has no domain instance")

    membership = conn.execute(
        "SELECT method, split, COUNT(*) FROM ddi_split_membership GROUP BY method, split"
    ).fetchall()
    got = {(m, s) for m, s, _n in membership}
    want = {(m, s) for m, splits in EXPECTED.items() for s in splits}
    check(got == want, f"ddi_split_membership covers {sorted(got)}, expected {sorted(want)}")
    for method, split, n in membership:
        check(n > 0, f"{method}/{split} has no membership rows")

    # Instance pairs of the external test split must have distinct parent proteins.
    bad = conn.execute(
        "SELECT COUNT(*) FROM ddi_split_membership AS m "
        "JOIN domain_protein_map AS a ON a.instance_id = m.instance_id_a "
        "JOIN domain_protein_map AS b ON b.instance_id = m.instance_id_b "
        "WHERE m.split = 'test' AND m.method LIKE 'external_test%' AND a.protein_id = b.protein_id"
    ).fetchone()[0]
    check(bad == 0, f"{bad} external-test instance pairs share a parent protein")

    external_test_ddis = conn.execute(
        "SELECT DISTINCT ddi.source FROM ddi_split_membership AS m "
        "JOIN domain_domain_interaction AS ddi ON ddi.id = m.ddi_id "
        "WHERE m.split = 'test' AND m.method LIKE 'external_test%'"
    ).fetchall()
    for (source,) in external_test_ddis:
        sources = sources_of(source)
        check("3did" not in sources and "sampled_negative" not in sources,
              f"external test split holds a protected-source DDI ({source})")
    if not external_test_ddis:
        notes.append("the external test split is empty: no external DDI survived the fixture's "
                     "family selection")
    conn.close()
    return {tuple(sorted((a, b))) for a, b, _n, _s in rows}


def check_splits(outdir, found, master_pairs):
    for method, splits in sorted(found.items()):
        for split in splits:
            path = os.path.join(outdir, "databases", method, f"{split}.sqlite3")
            conn = sqlite3.connect(path)
            rows = ddi_rows(conn)
            pairs = {tuple(sorted((a, b))) for a, b, _n, _s in rows}
            check(pairs <= master_pairs,
                  f"{method}/{split} holds {len(pairs - master_pairs)} DDIs absent from the master DB")
            in_split = conn.execute(
                "SELECT COUNT(*) FROM domain_domain_interaction AS ddi WHERE ddi.id NOT IN "
                "(SELECT ddi_id FROM ddi_split_membership WHERE method = ? AND split = ?)",
                (method, split),
            ).fetchone()[0]
            check(in_split == 0, f"{method}/{split} holds {in_split} DDIs not in that split")
            check(bool(rows), f"{method}/{split} is empty")
            foreign = conn.execute(
                "SELECT COUNT(*) FROM ddi_interaction_confirmed WHERE method != ?", (method,)
            ).fetchone()[0]
            check(foreign == 0, f"{method}/{split} carries {foreign} other methods' confirmation rows")
            check_types(conn, f"{method}/{split}")
            conn.close()


def check_types(conn, label):
    """Storage classes that consumers compare numerically."""
    bad_score = conn.execute(
        "SELECT COUNT(*) FROM protein_protein_interaction WHERE typeof(score) != 'real'"
    ).fetchone()[0]
    check(bad_score == 0, f"{label}: {bad_score} protein_protein_interaction.score values are not REAL")
    bad_pos = conn.execute(
        "SELECT COUNT(*) FROM domain_protein_map "
        "WHERE typeof(start_pos) != 'integer' OR typeof(end_pos) != 'integer'"
    ).fetchone()[0]
    check(bad_pos == 0, f"{label}: {bad_pos} domain_protein_map rows have non-INTEGER positions")
    bad_ddi = conn.execute(
        "SELECT COUNT(*) FROM ddi_interaction_confirmed WHERE typeof(ddi_id) != 'integer'"
    ).fetchone()[0]
    check(bad_ddi == 0, f"{label}: {bad_ddi} ddi_interaction_confirmed.ddi_id values are not INTEGER")


def check_family_disjointness(master):
    """minimal_leakage*: no Pfam family in two of a method's splits (F6). The two
    test splits are one test set under two negative sets, so they are one side."""
    conn = sqlite3.connect(master)
    rows = conn.execute("""
        WITH fam AS (
          SELECT DISTINCT m.method,
                 CASE WHEN m.split LIKE 'test%' THEN 'test' ELSE m.split END AS side,
                 d.pfam_id AS pfam
          FROM ddi_split_membership m
          JOIN domain_domain_interaction x ON x.id = m.ddi_id
          JOIN domain d ON d.id IN (x.domain_id_a, x.domain_id_b)
          WHERE m.method LIKE 'minimal_leakage%')
        SELECT a.method, a.side, b.side, COUNT(*)
        FROM fam a JOIN fam b ON a.method = b.method AND a.pfam = b.pfam AND a.side < b.side
        GROUP BY 1, 2, 3""").fetchall()
    conn.close()
    for method, side_a, side_b, n in rows:
        check(False, f"{method}: {n} Pfam families are in both {side_a} and {side_b}")


def check_manifest(outdir, master):
    path = os.path.join(outdir, "reports", "run_manifest.json")
    if not os.path.exists(path):
        return  # reported as a missing report in main()
    with open(path) as fh:
        manifest = json.load(fh).get("split_counts", {})
    conn = sqlite3.connect(master)
    for method, split, negative, n in conn.execute(
        "SELECT m.method, m.split, d.negative, COUNT(DISTINCT m.ddi_id) "
        "FROM ddi_split_membership m JOIN domain_domain_interaction d ON d.id = m.ddi_id "
        "GROUP BY 1, 2, 3"
    ):
        key = "negatives" if negative else "positives"
        got = manifest.get(method, {}).get(split, {}).get(key)
        check(got == n, f"run_manifest.json says {method}/{split} has {got} {key}, the master has {n}")
    conn.close()


def structure_key(instance_id_a, instance_id_b):
    """structures.h5 key -- mirrors utils_struct.structure_h5_keys()."""
    return (f"{instance_id_a.split('_')[0]}_{instance_id_b.split('_')[0]}",
            f"{instance_id_a}__{instance_id_b}")


def check_structure(outdir, master):
    conn = sqlite3.connect(master)
    cols = {r[1] for r in conn.execute("PRAGMA table_info(ddi_split_membership)")}
    check({"has_structure", "z_score"} <= cols,
          f"ddi_split_membership lacks has_structure/z_score (has {sorted(cols)})")
    if not {"has_structure", "z_score"} <= cols:
        conn.close()
        return

    bad_flag = conn.execute(
        "SELECT COUNT(*) FROM ddi_split_membership WHERE has_structure NOT IN (0, 1)"
    ).fetchone()[0]
    check(bad_flag == 0, f"{bad_flag} membership rows have has_structure outside 0/1")
    z_without = conn.execute(
        "SELECT COUNT(*) FROM ddi_split_membership WHERE z_score IS NOT NULL AND has_structure = 0"
    ).fetchone()[0]
    check(z_without == 0, f"{z_without} membership rows carry a z_score but no structure")
    # One score per (method, instance pair), copied to every split row of it.
    z_split = conn.execute(
        "SELECT COUNT(*) FROM (SELECT 1 FROM ddi_split_membership WHERE has_structure = 1 "
        "GROUP BY method, ddi_id, instance_id_a, instance_id_b "
        "HAVING COUNT(DISTINCT z_score) > 1 "
        "    OR (COUNT(z_score) > 0 AND COUNT(z_score) < COUNT(*)))"
    ).fetchone()[0]
    check(z_split == 0, f"{z_split} (method, instance pair)s have different z-scores across splits")
    # has_structure is a property of the instance pair, not of the method.
    flag_split = conn.execute(
        "SELECT COUNT(*) FROM (SELECT 1 FROM ddi_split_membership "
        "GROUP BY ddi_id, instance_id_a, instance_id_b HAVING COUNT(DISTINCT has_structure) > 1)"
    ).fetchone()[0]
    check(flag_split == 0, f"{flag_split} instance pairs have has_structure set on some rows only")

    unscored_votes = conn.execute(
        "SELECT COUNT(*) FROM ddi_interaction_confirmed c WHERE NOT EXISTS ("
        "  SELECT 1 FROM ddi_split_membership m WHERE m.ddi_id = c.ddi_id "
        "  AND m.method = c.method AND m.z_score IS NOT NULL)"
    ).fetchone()[0]
    check(unscored_votes == 0, f"{unscored_votes} ddi_interaction_confirmed rows have no scored instance pair")
    bad_vote = conn.execute(
        "SELECT COUNT(*) FROM ddi_interaction_confirmed "
        "WHERE majority_confirmed NOT IN (0, 1) OR mean_confirmed NOT IN (0, 1)"
    ).fetchone()[0]
    check(bad_vote == 0, f"{bad_vote} ddi_interaction_confirmed rows hold a non-boolean vote")

    triples = conn.execute(
        "SELECT DISTINCT instance_id_a, instance_id_b FROM ddi_split_membership WHERE has_structure = 1"
    ).fetchall()
    scored_methods = [m for (m,) in conn.execute("SELECT DISTINCT method FROM ddi_interaction_confirmed")]
    conn.close()

    h5_path = os.path.join(outdir, "structures", "structures.h5")
    if not triples:
        notes.append("no membership row has a structure: structure was off, or no "
                     "af_metadata pair matched -- structure checks skipped")
        return
    check(os.path.exists(h5_path), f"{len(triples)} instance pairs have a structure but {h5_path} is missing")
    try:
        import h5py
    except ImportError:
        notes.append("h5py is not installed, so structures.h5 was not opened")
        h5py = None
    if h5py is not None and os.path.exists(h5_path):
        with h5py.File(h5_path, "r") as h5:
            want = {structure_key(a, b) for a, b in triples}
            have = {(g, d) for g in h5 for d in h5[g]}
            check(not (want - have), f"{len(want - have)} has_structure = 1 pairs are missing from structures.h5")
            check(not (have - want), f"structures.h5 holds {len(have - want)} slices no membership row claims")
            for group, dset in sorted(want & have)[:20]:
                raw = h5[group][dset][()].tobytes()
                check(raw[:2] == b"\x1f\x8b", f"structures.h5 {group}/{dset} is not gzip data")

    for method in scored_methods:
        check_scoring_matrix(outdir, method)


def check_scoring_matrix(outdir, method):
    root = os.path.join(outdir, "ddi_scoring")
    paths = {k: os.path.join(root, f"{method}.{k}") for k in ("c_ab_matrix.csv", "db_freq.csv", "t_db.txt")}
    missing = [k for k, v in paths.items() if not os.path.exists(v)]
    check(not missing, f"ddi_scoring/{method}: missing {missing}")
    if missing:
        return
    with open(paths["c_ab_matrix.csv"]) as fh:
        c_ab = {(r["aaA"], r["aaB"]): int(r["count"]) for r in csv.DictReader(fh)}
    check(len(c_ab) == 400, f"ddi_scoring/{method}.c_ab_matrix.csv has {len(c_ab)} cells, expected 400")
    asym = [k for k, v in c_ab.items() if c_ab.get((k[1], k[0])) != v]
    check(not asym, f"ddi_scoring/{method}.c_ab_matrix.csv is not symmetric ({len(asym)} cells)")
    with open(paths["db_freq.csv"]) as fh:
        total = sum(float(r["frequency"]) for r in csv.DictReader(fh))
    check(abs(total - 1.0) < 1e-6 or total == 0,
          f"ddi_scoring/{method}.db_freq.csv sums to {total}")
    with open(paths["t_db.txt"]) as fh:
        t_db = int(fh.read().strip())
    aas = sorted({a for a, _ in c_ab})
    expected = sum(c_ab[(a, a)] for a in aas) + 2 * sum(
        c_ab[(a, b)] for i, a in enumerate(aas) for b in aas[i + 1:])
    check(t_db == expected, f"ddi_scoring/{method}.t_db.txt is {t_db}, the matrix gives {expected}")


def check_embeddings(outdir, master):
    """Check every published embedding file against the benchmark's key contract."""
    root = os.path.join(outdir, "embeddings")
    if not os.path.isdir(root):
        notes.append("no embeddings/ directory: expected only when "
                     "--embedding_models_domainbench is empty")
        return
    files = sorted(f for f in os.listdir(root) if f.endswith(".h5"))
    check(bool(files), "embeddings/ exists but holds no .h5 file")
    if not files:
        return
    try:
        import h5py
    except ImportError:
        notes.append("h5py is not installed, so the embedding files were not opened")
        return

    conn = sqlite3.connect(master)
    instances = conn.execute(
        "SELECT d.pfam_id, COALESCE(m.instance_id, 'r' || m.rowid) "
        "FROM domain_protein_map AS m JOIN domain AS d ON d.id = m.domain_id"
    ).fetchall()
    conn.close()
    check(bool(instances), "the master database has no domain instances at all")

    for name in files:
        path = os.path.join(root, name)
        with h5py.File(path, "r") as h5:
            for attr in EMBEDDING_ATTRS:
                check(attr in h5.attrs, f"{name} is missing the root attribute {attr!r}")
            # Not every instance has to resolve -- --embedding_max_len drops the
            # long tail -- but the *shape* has to, or the benchmark skips silently.
            resolved = sum(
                1 for pfam_id, key in instances
                if pfam_id in h5 and key in h5[pfam_id]
            )
            check(resolved > 0,
                  f"{name} resolves none of the master's {len(instances)} "
                  "(pfam_id, instance_key) pairs -- the key layout has drifted")
            if 0 < resolved < len(instances):
                notes.append(f"{name} covers {resolved} of {len(instances)} domain instances "
                             "(the rest were dropped by --embedding_max_len or an OOM)")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    master = os.path.join(args.outdir, "domainsplit.sqlite3")
    check(os.path.exists(master), f"{master} does not exist")
    for report in ("reports/source_conflicts.tsv", "reports/pruned_ddis.tsv",
                   "reports/external_test_dropped.tsv",
                   "reports/ddi_source_attrition.tsv",
                   "reports/ddi_tier_breakdown.tsv",
                   "reports/run_manifest.json"):
        check(os.path.exists(os.path.join(args.outdir, report)), f"missing {report}")

    found = check_layout(args.outdir)
    if os.path.exists(master):
        master_pairs = check_master(master)
        check_splits(args.outdir, found, master_pairs)
        conn = sqlite3.connect(master)
        check_types(conn, "master")
        conn.close()
        check_family_disjointness(master)
        check_manifest(args.outdir, master)
        check_structure(args.outdir, master)
        check_embeddings(args.outdir, master)

    for note in notes:
        print(f"[check] note: {note}")
    if failures:
        for failure in failures:
            print(f"[check] FAIL: {failure}", file=sys.stderr)
        sys.exit(1)
    print(f"[check] ok: {sum(len(s) for s in found.values())} split databases")


if __name__ == "__main__":
    main()
