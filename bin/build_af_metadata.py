#!/usr/bin/env python3
"""Build af_metadata.csv, the index of AF complex models ENRICH_STRUCTURAL_AF
reads, from three optional sources:

  --manifest  the raw AFDB-heterodimer download manifest
              (uniprot_id_a,uniprot_id_b,model_entity_id,local_tar_name,path;
              `path` is a directory holding <model_entity_id>-model_v1.pdb.zst),
              source `downloaded`
  --inferred_dir  a root scanned recursively for AF3 job directories
              <name>/<name>_model.cif with <name> = <prefix>_<uniprot_a>_<uniprot_b>
              (the job builder's default prefix is `Pos_` for every pair, so the
              prefix is never read as a label), source `inferred`; repeatable
  --af_metadata   an existing CSV uniprot_id_a,uniprot_id_b,path,source (a
              relative path is relative to the CSV's real location)

Written from scratch on every run, never appended. Rules, in this order:
  1. A candidate whose model file does not exist is dropped (status model_missing).
  2. Exact repeats of one model (same unordered pair, source and model path) are
     counted once.
  3. Competing downloaded models are ranked by mean C-alpha pLDDT. Ties for the
     highest score and conflicting inferred models are fatal.
  4. A pair with a usable `downloaded` model keeps it; its `inferred` model is
     `shadowed`. Validation comes before precedence, so a missing downloaded file
     never hides a usable inferred one.

Each selected row carries a stat fingerprint (size, mtime_ns) of the model and of
its confidence sidecars (-1 for a missing sidecar). The model files are read from
paths inside this CSV and are not staged by Nextflow, so the fingerprint is what
makes a replaced model or sidecar change the CSV and, through its digest,
re-run ENRICH_STRUCTURAL_AF. It is a stat fingerprint, not a content hash.
The workflow runs this twice: an uncached scan whose --digest keys the cached
build (--expect_digest), so an unchanged scan reuses the build and everything
downstream of it.

Outputs: --out (the CSV), --inventory (every candidate with its status),
--report (counts per source), --digest (SHA-256 of --out).
"""
import argparse
import csv
import hashlib
import os
import re
import sys
from statistics import mean

from enrich_structural_af import confidence_paths, load_structure

OUT_COLS = ["uniprot_id_a", "uniprot_id_b", "path", "source",
            "model_size", "model_mtime_ns", "pae_size", "pae_mtime_ns",
            "summary_size", "summary_mtime_ns"]
# UniProt accession (optionally an isoform, -<n>): the job name's two last tokens
# must be accessions, so `batch_2026_11` is not mistaken for a job.
ACC = r"(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-[0-9]+)?"
JOB_DIR_RE = re.compile(rf"^(?P<prefix>[A-Za-z0-9]+)_(?P<a>{ACC})_(?P<b>{ACC})$")
SOURCES = ("downloaded", "inferred")


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifest")
    p.add_argument("--inferred_dir", action="append", default=[])
    p.add_argument("--af_metadata")
    p.add_argument("--out", required=True)
    p.add_argument("--inventory", required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--digest", required=True)
    p.add_argument("--expect_digest",
                   help="fail unless --out hashes to this (the uncached scan's digest): a mismatch "
                        "means the inferred directories changed between scan and build")
    return p.parse_args(argv)


def stat(path):
    """(size, mtime_ns), or (-1, -1) if absent / not derivable."""
    if not path:
        return -1, -1
    try:
        st = os.stat(path)
    except FileNotFoundError:
        return -1, -1
    return st.st_size, st.st_mtime_ns


def read_manifest(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    need = {"uniprot_id_a", "uniprot_id_b", "model_entity_id", "path"}
    if rows and not need <= set(rows[0]):
        sys.exit(f"--manifest {path}: needs columns {sorted(need)}, has {sorted(rows[0])}")
    for i, r in enumerate(rows, start=2):
        if not all((r.get(c) or "").strip() for c in need):
            sys.exit(f"--manifest {path}: empty field on row {i}")
        yield (r["uniprot_id_a"].strip(), r["uniprot_id_b"].strip(),
               os.path.join(os.path.abspath(r["path"].strip()), f"{r['model_entity_id'].strip()}-model_v1.pdb.zst"),
               "downloaded", f"manifest:{i}")


def scan_inferred(root):
    if not os.path.isdir(root):
        sys.exit(f"--inferred_dir {root}: not a directory")
    root = os.path.abspath(root)
    m = JOB_DIR_RE.match(os.path.basename(root))
    if m:  # the root is itself one job directory
        name = os.path.basename(root)
        yield (m["a"], m["b"], os.path.join(root, f"{name}_model.cif"), "inferred", f"{root}:.")
        return
    for dirpath, dirnames, _files in os.walk(root):
        jobs = [name for name in dirnames if JOB_DIR_RE.match(name)]
        # A job directory holds per-seed subdirectories, never another job: do
        # not descend (tens of thousands of directories on NFS).
        dirnames[:] = sorted(set(dirnames) - set(jobs))
        for name in sorted(jobs):
            m = JOB_DIR_RE.match(name)
            yield (m["a"], m["b"], os.path.join(dirpath, name, f"{name}_model.cif"),
                   "inferred", f"{root}:{os.path.relpath(os.path.join(dirpath, name), root)}")


def read_af_metadata(path):
    base = os.path.dirname(os.path.realpath(path))
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    need = {"uniprot_id_a", "uniprot_id_b", "path", "source"}
    if rows and not need <= set(rows[0]):
        sys.exit(f"--af_metadata {path}: needs columns {sorted(need)}, has {sorted(rows[0])}")
    for i, r in enumerate(rows, start=2):
        if not all((r.get(c) or "").strip() for c in need):
            sys.exit(f"--af_metadata {path}: empty field on row {i}")
        if r["source"] not in SOURCES:
            sys.exit(f"--af_metadata {path}: row {i} has source {r['source']!r}, expected one of {SOURCES}")
        yield (r["uniprot_id_a"].strip(), r["uniprot_id_b"].strip(),
               os.path.join(base, r["path"].strip()), r["source"], f"af_metadata:{i}")


def mean_plddt(path):
    structure, _ = load_structure(path, "downloaded")
    scores = [residue["CA"].bfactor for residue in structure[0].get_residues() if "CA" in residue]
    if not scores or not all(0 <= score <= 100 for score in scores):
        sys.exit(f"{path}: missing or invalid C-alpha pLDDT")
    return mean(scores)


def build(candidates):
    """-> (selected rows, inventory rows, report rows). Exits on a conflict."""
    inventory, usable = [], {}   # usable[(key, source)] = {path: (a, b, origin)}
    counts = {s: {"candidates": 0, "model_missing": 0, "exact_duplicate": 0,
                  "shadowed": 0, "selected": 0} for s in SOURCES}
    for a, b, path, source, origin in candidates:
        counts[source]["candidates"] += 1
        key = tuple(sorted((a, b)))
        if not os.path.isfile(path):
            counts[source]["model_missing"] += 1
            inventory.append([a, b, source, path, "model_missing", origin])
            continue
        seen = usable.setdefault((key, source), {})
        if path in seen:
            counts[source]["exact_duplicate"] += 1
            inventory.append([a, b, source, path, "exact_duplicate", origin])
            continue
        seen[path] = (a, b, origin)

    conflicts = []
    for (key, source), paths in usable.items():
        if len(paths) == 1:
            continue
        if source == "downloaded":
            ranked = sorted(((mean_plddt(path), path) for path in paths), reverse=True)
            if ranked[0][0] > ranked[1][0]:
                for _, path in ranked[1:]:
                    a, b, origin = paths.pop(path)
                    counts[source]["shadowed"] += 1
                    inventory.append([a, b, source, path, "shadowed", origin])
                continue
        conflicts.append((key, source, sorted(paths)))
    if conflicts:
        lines = "\n".join(f"  {'/'.join(k)} ({s}): {', '.join(p)}" for k, s, p in conflicts[:10])
        sys.exit(f"{len(conflicts)} protein pair(s) have two different {conflicts[0][1]} models; "
                 f"no rule says which is right, remove one:\n{lines}")

    selected = []
    for (key, source), paths in sorted(usable.items()):
        (path, (a, b, origin)), = paths.items()
        if source == "inferred" and (key, "downloaded") in usable:
            counts[source]["shadowed"] += 1
            inventory.append([a, b, source, path, "shadowed", origin])
            continue
        counts[source]["selected"] += 1
        inventory.append([a, b, source, path, "selected", origin])
        summary, pae = confidence_paths(path, source)
        selected.append([a, b, path, source, *stat(path), *stat(pae), *stat(summary)])
    report = [[s, *counts[s].values()] for s in SOURCES]
    return selected, sorted(inventory), report


def write_tsv(path, header, rows, sep="\t"):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter=sep, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def main(argv=None):
    args = parse_args(argv)
    if not (args.manifest or args.inferred_dir or args.af_metadata):
        sys.exit("give at least one of --manifest, --inferred_dir, --af_metadata")

    def candidates():
        if args.af_metadata:
            yield from read_af_metadata(args.af_metadata)
        if args.manifest:
            yield from read_manifest(args.manifest)
        for root in args.inferred_dir:
            yield from scan_inferred(root)

    selected, inventory, report = build(candidates())
    if not selected:
        sys.exit("no usable AF model found in the given sources")
    write_tsv(args.out, OUT_COLS, selected, sep=",")
    write_tsv(args.inventory, ["uniprot_id_a", "uniprot_id_b", "source", "path", "status", "origin"], inventory)
    write_tsv(args.report, ["source", "candidates", "model_missing", "exact_duplicate", "shadowed", "selected"], report)
    with open(args.out, "rb") as fh:
        digest = hashlib.sha256(fh.read()).hexdigest()
    with open(args.digest, "w") as fh:
        fh.write(digest + "\n")
    if args.expect_digest and args.expect_digest.strip() != digest:
        sys.exit(f"af_metadata changed between the scan ({args.expect_digest.strip()}) and the build "
                 f"({digest}): AF jobs were added or replaced during this run. Re-run with -resume.")
    for row in report:
        print("[build_af_metadata] " + " ".join(f"{h}={v}" for h, v in zip(
            ["source", "candidates", "model_missing", "exact_duplicate", "shadowed", "selected"], row)))


if __name__ == "__main__":
    main()
