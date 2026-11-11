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
* ``--db``: the master; per-(method, split) positive/negative counts are read
  from ``ddi_split_membership``, so the class ratio every published split
  *actually* has is on record, not the ratio it was asked for.

URLs get a HEAD request for ``Last-Modified``/``ETag``. That is the server's
state when the manifest is written, which is the download's state unless a file
came from ``cache_dir`` -- the reason ``cache_dir`` is recorded alongside.
A failed HEAD is recorded as an error, never fatal: provenance must not be able
to fail a finished run.
"""

import argparse
import hashlib
import json
import sqlite3
import sys
import urllib.request


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
            "SELECT m.method, m.split, d.negative, COUNT(DISTINCT m.ddi_id) "
            "FROM ddi_split_membership m JOIN domain_domain_interaction d ON d.id = m.ddi_id "
            "GROUP BY 1, 2, 3"
        ).fetchall()
    finally:
        conn.close()
    out = {}
    for method, split, negative, n in rows:
        c = out.setdefault(method, {}).setdefault(split, {"positives": 0, "negatives": 0})
        c["negatives" if negative else "positives"] += n
    for splits in out.values():
        for c in splits.values():
            c["negatives_per_positive"] = round(c["negatives"] / c["positives"], 4) if c["positives"] else None
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--base", required=True)
    p.add_argument("--input", action="append", default=[], metavar="NAME=PATH")
    p.add_argument("--db", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--no-network", action="store_true", help="skip the HEAD requests (tests, offline runs)")
    p.add_argument("--versions")
    p.add_argument("--process-name", default="WRITE_RUN_MANIFEST")
    args = p.parse_args(argv)

    with open(args.base) as f:
        manifest = json.load(f)

    inputs = {}
    for spec in args.input:
        name, _, path = spec.partition("=")
        if not path:
            sys.exit(f"--input expects NAME=PATH, got {spec!r}")
        inputs[name] = {"path": path, "sha256": sha256(path)}
    manifest["manual_inputs"] = inputs

    urls = manifest.get("urls", {})
    manifest["url_headers"] = {
        name: ({"url": u} if args.no_network else {"url": u, **head(u)})
        for name, u in urls.items()
    }
    manifest["split_counts"] = split_counts(args.db)
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
