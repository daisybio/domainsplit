#!/usr/bin/env python3
"""Write ``reports/run_manifest.json``: what a run was made from and what it made.

The 2026-09 audit (DATASET_AUDIT.md F10) could not establish what produced the
published splits: the work directories were gone, the three manual inputs are
unversioned paths, and every ``url_*`` points at a "current" release. This file
is the record that makes the next audit answerable from ``results/`` alone.

Inputs:

* ``--base``: JSON the workflow builds in Groovy -- every param, the git
  commits of this repo and the ppi-splitting submodule, Nextflow/session info,
  and the URLs the run downloads from.
* ``--input NAME=PATH`` (repeatable): manual inputs, hashed with SHA-256.
* ``--local NAME=PATH`` (repeatable): a ``urls`` entry given as a local path,
  staged by the workflow; hashed into ``url_headers[NAME]`` instead of HEADed.
* ``--db``: the master; per-(method, split) positive/negative counts are read
  from ``ddi_split_membership``, so the class ratio every published split
  *actually* has is on record, not the ratio it was asked for. Two units:
  ``positives``/``negatives``/``negatives_per_positive`` count distinct DDIs
  (Pfam family pairs, what ppi-splitting balances), while
  ``positive_instance_pairs``/``negative_instance_pairs``/
  ``instance_pairs_negatives_per_positive`` count membership rows (the concrete
  examples a split DB holds); they differ when negative expansion falls short.
  ``sampled_negatives_without_membership`` counts protected ``sampled_negative``
  DDIs in the master with no membership row: they are in no split, yet still
  block an external label for their pair.

URLs get a HEAD request for ``Last-Modified``/``ETag``. That is the server's
state when the manifest is written, which is the download's state unless a file
came from ``cache_dir`` -- the reason ``cache_dir`` is recorded alongside.
A failed HEAD is recorded as an error, never fatal: provenance must not be able
to fail a finished run.
"""

import argparse
import hashlib
import json
import os
import sqlite3
import sys
import urllib.request


# BUILD_AF3_QUEUE / STRUCTURE_COVERAGE_BALANCE / FETCH_INTERPRO_METADATA
# summaries arriving via --input.
EMBEDDED = {"af3_queue.json": "af3_queue",
            "structure_coverage_balance.json": "structure_coverage_balance",
            "interpro_metadata.json": "interpro_metadata"}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def head(url, timeout=30):
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return {k: r.headers.get(k) for k in ("Last-Modified", "ETag", "Content-Length")}
    except Exception as e:  # noqa: BLE001 -- recorded, never fatal
        return {"error": f"{type(e).__name__}: {e}"}


def split_counts(db):
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        rows = conn.execute(
            "SELECT m.method, m.split, d.negative, COUNT(DISTINCT m.ddi_id), COUNT(*) "
            "FROM ddi_split_membership m JOIN domain_domain_interaction d ON d.id = m.ddi_id "
            "GROUP BY 1, 2, 3"
        ).fetchall()
    finally:
        conn.close()
    out = {}
    for method, split, negative, n, n_pairs in rows:
        c = out.setdefault(method, {}).setdefault(
            split, {"positives": 0, "negatives": 0, "positive_instance_pairs": 0, "negative_instance_pairs": 0}
        )
        c["negatives" if negative else "positives"] += n
        c["negative_instance_pairs" if negative else "positive_instance_pairs"] += n_pairs
    for splits in out.values():
        for c in splits.values():
            c["negatives_per_positive"] = round(c["negatives"] / c["positives"], 4) if c["positives"] else None
            c["instance_pairs_negatives_per_positive"] = (
                round(c["negative_instance_pairs"] / c["positive_instance_pairs"], 4)
                if c["positive_instance_pairs"] else None
            )
    return out


def membership_fingerprint(db):
    """SHA-256 over every membership row as (method, split, pfam_a, pfam_b,
    instance_a, instance_b, negative), sorted -- no surrogate ids, so two runs
    with the same splits get the same value. Equal fingerprint = AF3 queues and
    predictions made against one run are valid for the other."""
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    h, n = hashlib.sha256(), 0
    try:
        for row in conn.execute(
            "SELECT m.method, m.split, da.pfam_id, db.pfam_id, COALESCE(m.instance_id_a, ''), "
            "COALESCE(m.instance_id_b, ''), CAST(d.negative AS INTEGER) "
            "FROM ddi_split_membership m JOIN domain_domain_interaction d ON d.id = m.ddi_id "
            "JOIN domain da ON da.id = d.domain_id_a JOIN domain db ON db.id = d.domain_id_b "
            "ORDER BY 1, 2, 3, 4, 5, 6, 7"
        ):
            h.update(("\t".join(map(str, row)) + "\n").encode())
            n += 1
    finally:
        conn.close()
    return {"sha256": h.hexdigest(), "rows": n}


def protein_fingerprint(db):
    """SHA-256 over every (uniprot_id, sequence) of `protein`, sorted: a frozen
    run reuses the splits' domain coordinates, which are only valid on the
    same parent sequences."""
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    h, n = hashlib.sha256(), 0
    try:
        for uniprot_id, sequence in conn.execute(
                "SELECT uniprot_id, COALESCE(sequence, '') FROM protein ORDER BY uniprot_id"):
            h.update(f"{uniprot_id}\t{sequence}\n".encode())
            n += 1
    finally:
        conn.close()
    return {"sha256": h.hexdigest(), "rows": n}


def sampled_negatives_without_membership(db):
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        return conn.execute(
            "SELECT COUNT(*) FROM domain_domain_interaction d "
            "WHERE d.negative = 1 "
            "AND INSTR(',' || REPLACE(d.source, ' ', '') || ',', ',sampled_negative,') > 0 "
            "AND NOT EXISTS (SELECT 1 FROM ddi_split_membership m WHERE m.ddi_id = d.id)"
        ).fetchone()[0]
    finally:
        conn.close()


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--base", required=True)
    p.add_argument("--input", action="append", default=[], metavar="NAME=PATH")
    p.add_argument("--local", action="append", default=[], metavar="NAME=PATH",
                   help="a urls entry that is a local file: hashed, not HEADed (also with --no-network)")
    p.add_argument("--db", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--no-network", action="store_true", help="skip the HEAD requests (tests, offline runs)")
    p.add_argument("--versions")
    p.add_argument("--process-name", default="WRITE_RUN_MANIFEST")
    args = p.parse_args(argv)

    with open(args.base) as f:
        manifest = json.load(f)

    def pairs(flag, specs):
        out = {}
        for spec in specs:
            name, _, path = spec.partition("=")
            if not path:
                sys.exit(f"{flag} expects NAME=PATH, got {spec!r}")
            # A staged file's name can carry its stageAs directory (`extra/x.tsv`).
            out[os.path.basename(name)] = path
        return out

    manifest["manual_inputs"] = {
        name: {"path": path, "sha256": sha256(path)} for name, path in pairs("--input", args.input).items()
    }

    # Step summaries hashed above and also embedded under their manifest key.
    for name, key in EMBEDDED.items():
        path = pairs("--input", args.input).get(name)
        if path:
            with open(path) as f:
                manifest[key] = json.load(f)

    local = pairs("--local", args.local)
    urls = manifest.get("urls", {})
    manifest["url_headers"] = {
        name: ({"url": u, "sha256": sha256(local[name])} if name in local
               else {"url": u} if args.no_network else {"url": u, **head(u)})
        for name, u in urls.items()
    }
    manifest["split_counts"] = split_counts(args.db)
    manifest["sampled_negatives_without_membership"] = sampled_negatives_without_membership(args.db)
    manifest["membership_fingerprint"] = membership_fingerprint(args.db)
    manifest["protein_fingerprint"] = protein_fingerprint(args.db)
    # ppi-splitting publishes no solver status, so whether an ILP hit its time
    # limit (and returned an incumbent rather than an optimum) is not knowable here.
    manifest["ilp_solver_status"] = "not recorded: ppi-splitting does not publish it"

    with open(args.out, "w") as f:
        json.dump(manifest, f, indent=2, sort_keys=True)
        f.write("\n")

    if args.versions:
        with open(args.versions, "w") as f:
            f.write(f'"{args.process_name}":\n')
            f.write(f"    python: {sys.version.split()[0]}\n")
            f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
