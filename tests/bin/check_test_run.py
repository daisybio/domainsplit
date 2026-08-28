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
7. the reports exist;
8. every published embedding HDF5 satisfies the benchmark's own key contract --
   `h5[str(domain_id)][COALESCE(instance_id, 'r' || rowid)]` resolves for the
   master's domain instances -- and carries its root attributes. This is the seam
   that yields zero training rows *silently* if it drifts, so it is checked here
   as well as in tests/python/test_export_domain_embeddings.py.

Usage: ``tests/bin/check_test_run.py --outdir results``
"""

import argparse
import os
import sqlite3
import sys

# `params.ppi_splitting_multi_negset = true` (the default) asks ppi-splitting for
# both ILP negative sets per ILP row, so each of those rows yields a plain
# directory and an `_hcni` one: 5 method directories, 18 split databases. The
# `random` row stays `uniform`-only -- the naive baseline stays a baseline.
# Set `--ppi_splitting_multi_negset false` and this drops back to 3 / 11.
INTERNAL_SPLITS = ["train", "validation", "test_balanced", "test_realistic"]
EXTERNAL_SPLITS = ["train", "validation", "test"]
EXPECTED = {
    "random": INTERNAL_SPLITS,
    "minimal_leakage": INTERNAL_SPLITS,
    "minimal_leakage_hcni": INTERNAL_SPLITS,
    "external_test": EXTERNAL_SPLITS,
    "external_test_hcni": EXTERNAL_SPLITS,
}
EXTERNAL_SOURCES = ("single_domain_ppi", "PPIDM", "negatome")
# Root attributes export_domain_embeddings.py writes; `domainsplit_run` is the
# cross-run guard, since `domain.id` is a surrogate integer that means nothing
# outside the run that assigned it.
EMBEDDING_ATTRS = ("model", "pooling", "dim", "dtype", "key_layout",
                   "n_domains", "n_instances", "domainsplit_run")

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
            conn.close()


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
        "SELECT domain_id, COALESCE(instance_id, 'r' || rowid) FROM domain_protein_map"
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
                1 for domain_id, key in instances
                if str(domain_id) in h5 and key in h5[str(domain_id)]
            )
            check(resolved > 0,
                  f"{name} resolves none of the master's {len(instances)} "
                  "(domain_id, instance_key) pairs -- the key layout has drifted")
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
                   "reports/external_test_dropped.tsv"):
        check(os.path.exists(os.path.join(args.outdir, report)), f"missing {report}")

    found = check_layout(args.outdir)
    if os.path.exists(master):
        master_pairs = check_master(master)
        check_splits(args.outdir, found, master_pairs)
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
