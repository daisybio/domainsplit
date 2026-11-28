#!/usr/bin/env python3
"""Freeze a run's splits, and check a frozen bundle before it is reused.

`export` copies what a run handed from ppi-splitting to INGEST_SPLITS -- the
family- and instance-level split CSVs per (method, split) -- together with
FETCH_DOMAIN_META's outputs (the sampled domain instances and their sequences,
species, GO and dropped families) and the positive population EXPORT_SPLIT_DDIS
handed to ppi-splitting, into one directory with a bundle_manifest.json (SHA-256
of every file, the split list, the params the bundle depends on).

`verify` is the gate of a `--frozen_splits` run: every bundle file must still
hash to its manifest entry, and the positive population this run rebuilt from
its own (pinned) sources must equal the frozen one line for line, ignoring row
order. Anything else would ingest splits drawn from a different population.

`check_db` is the second gate, before structure and split databases: the
frozen run's final master must carry the same membership rows (external test
included) and the same protein sequences as the exporting run (fingerprints of
write_run_manifest.py, recorded by `export --db`). It catches what the first
gate cannot see: changed external-test inputs or params, split rows dropped at
ingest, a UniProt sequence changed under the same accession.

A frozen run reuses the splits; it does not reproduce how they were drawn.
"""
import argparse
import hashlib
import json
import os
import shutil
import sys

from write_run_manifest import membership_fingerprint, protein_fingerprint

FORMAT = 2  # 2: fingerprints of the final master are required


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_lines(path):
    """(header, sorted body lines): EXPORT_SPLIT_DDIS's row order is not a contract."""
    with open(path) as fh:
        lines = fh.read().splitlines()
    if not lines:
        return "", []
    return lines[0], sorted(line for line in lines[1:] if line.strip())


def canonical_sha(path):
    header, body = canonical_lines(path)
    return hashlib.sha256("\n".join([header, *body]).encode()).hexdigest()


def spec(value, n):
    parts = value.split(":", n - 1)
    if len(parts) != n or not all(parts):
        sys.exit(f"expected {n} ':'-separated fields, got {value!r}")
    return parts


def export(args):
    out = args.out
    if os.path.exists(out):
        sys.exit(f"{out} exists")
    files, manifest = {}, {"format": FORMAT, "fetch": {}, "family": [], "instance": [],
                           "params": json.loads(args.params)}

    def put(src, rel):
        dst = os.path.join(out, rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copyfile(src, dst)
        files[rel] = sha256(dst)
        return rel

    for path in args.fetch:
        name = os.path.basename(path)
        if name in manifest["fetch"]:
            sys.exit(f"--fetch: two files named {name}")
        manifest["fetch"][name] = put(path, f"fetch/{name}")
    manifest["split_ddis"] = put(args.split_ddis, f"split_ddis/{os.path.basename(args.split_ddis)}")
    manifest["split_ddis_canonical_sha256"] = canonical_sha(args.split_ddis)
    for level in ("family", "instance"):
        seen = set()
        for value in getattr(args, level):
            method, split, path = spec(value, 3)
            if (method, split) in seen:
                sys.exit(f"--{level}: {method}/{split} given twice")
            seen.add((method, split))
            rel = put(path, f"splits/{level}/{method}/{split}.csv")
            manifest[level].append({"method": method, "split": split, "path": rel})
        manifest[level].sort(key=lambda e: (e["method"], e["split"]))
    if {(e["method"], e["split"]) for e in manifest["family"]} != \
            {(e["method"], e["split"]) for e in manifest["instance"]}:
        sys.exit("family and instance splits cover different (method, split) sets")
    manifest["files"] = dict(sorted(files.items()))
    manifest["fingerprints"] = {"membership": membership_fingerprint(args.db),
                                "protein": protein_fingerprint(args.db)}
    with open(os.path.join(out, "bundle_manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")


def verify(args):
    with open(os.path.join(args.bundle, "bundle_manifest.json")) as fh:
        manifest = json.load(fh)
    if manifest.get("format") != FORMAT:
        sys.exit(f"bundle format {manifest.get('format')!r}, this pipeline reads {FORMAT}")
    referenced = ({manifest["split_ddis"]} | set(manifest["fetch"].values())
                  | {e["path"] for level in ("family", "instance") for e in manifest[level]})
    unlisted = sorted(referenced - set(manifest["files"]))
    if unlisted:
        sys.exit(f"bundle references file(s) without a hash: {unlisted[:10]}")
    if {(e["method"], e["split"]) for e in manifest["family"]} != \
            {(e["method"], e["split"]) for e in manifest["instance"]}:
        sys.exit("bundle family and instance splits cover different (method, split) sets")
    bad = [rel for rel, digest in manifest["files"].items()
           if not os.path.isfile(os.path.join(args.bundle, rel))
           or sha256(os.path.join(args.bundle, rel)) != digest]
    if bad:
        sys.exit(f"{len(bad)} bundle file(s) missing or changed since export: {bad[:10]}")
    if "fingerprints" not in manifest:
        sys.exit("bundle has no master fingerprints; re-export it with this pipeline version")
    if not args.split_ddis:
        print(f"[frozen_splits] bundle integrity verified: {len(manifest['files'])} files")
        return
    if canonical_sha(args.split_ddis) != manifest["split_ddis_canonical_sha256"]:
        fh_, frozen = canonical_lines(os.path.join(args.bundle, manifest["split_ddis"]))
        h_, now = canonical_lines(args.split_ddis)
        only_frozen, only_now = sorted(set(frozen) - set(now)), sorted(set(now) - set(frozen))
        sys.exit("this run's positive population (EXPORT_SPLIT_DDIS) differs from the frozen one: "
                 f"header {'same' if fh_ == h_ else f'{fh_!r} vs {h_!r}'}, "
                 f"{len(only_frozen)} row(s) only frozen (e.g. {only_frozen[:3]}), "
                 f"{len(only_now)} only in this run (e.g. {only_now[:3]}). "
                 "Its sources or params changed; pin them to the frozen run's or drop --frozen_splits.")
    print(f"[frozen_splits] bundle verified: {len(manifest['files'])} files, "
          f"{len(manifest['family'])} (method, split) pairs")


def check_db(args):
    with open(os.path.join(args.bundle, "bundle_manifest.json")) as fh:
        frozen = json.load(fh)["fingerprints"]
    now = {"membership": membership_fingerprint(args.db), "protein": protein_fingerprint(args.db)}
    diff = [k for k in ("membership", "protein") if now[k] != frozen[k]]
    if diff:
        sys.exit("this run's master differs from the frozen run's: " + "; ".join(
            f"{k} rows {frozen[k]['rows']} -> {now[k]['rows']}, sha256 {frozen[k]['sha256'][:12]} -> {now[k]['sha256'][:12]}"
            for k in diff) + ". AF3 queues and predictions made for the frozen run do not carry over; "
            "pin the sources (and external-test params) to the frozen run's.")
    print(f"[frozen_splits] master matches the frozen run: {now['membership']['rows']} membership rows, "
          f"{now['protein']['rows']} proteins")


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    e = sub.add_parser("export")
    e.add_argument("--fetch", action="append", default=[], metavar="PATH")
    e.add_argument("--db", required=True, help="the run's final master (fingerprints for check_db)")
    e.add_argument("--split_ddis", required=True)
    e.add_argument("--family", action="append", default=[], metavar="METHOD:SPLIT:PATH")
    e.add_argument("--instance", action="append", default=[], metavar="METHOD:SPLIT:PATH")
    e.add_argument("--params", required=True, help="JSON object of the params the splits depend on")
    e.add_argument("--out", required=True)
    v = sub.add_parser("verify")
    v.add_argument("--bundle", required=True)
    v.add_argument("--split_ddis", help="this run's positive population; omit to check bundle integrity only")
    c = sub.add_parser("check_db")
    c.add_argument("--bundle", required=True)
    c.add_argument("--db", required=True)
    args = p.parse_args(argv)
    {"export": export, "verify": verify, "check_db": check_db}[args.cmd](args)


if __name__ == "__main__":
    main()
