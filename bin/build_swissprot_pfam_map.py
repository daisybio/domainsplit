#!/usr/bin/env python3
"""Build a UniProt -> Pfam map for single-domain detection.

Reads the ``Entry/Entry Name/Gene Names/Pfam`` TSV that ``PARSE_SWISSPROT``
carves out of the SwissProt flat file (a URL is still accepted, for the same TSV
served over http) and emits ``swissprot_pfam_map.json``:

    {"accession_to_pfams": {accession: [Pfam, ...]}}

Only accessions. There used to be a ``name_to_accession`` map so the single-domain
step could resolve HIPPIE identifiers that were entry names or gene names, with
ambiguous gene names dropped; ``HIPPIE-current.txt`` carries UniProt accessions in
columns 1 and 4, so nothing needs it and an ambiguity that had to be resolved by
dropping rows is gone with it.

The Pfam lists come from UniProt's own ``DR Pfam`` cross-references rather than
from ``Pfam-A.regions.tsv.gz``, which is the other table that could answer "which
families does this protein carry". Those two agree closely -- both derive from the
same Pfam matches -- and reading regions here would mean a second pass over a
4.7 GB file at a point in the DAG that runs *before* the Pfam fetch, for a
disagreement that switching the instance source already removed.
"""

import argparse
import gzip
import json
import os
import shutil
import ssl
import sys
import urllib.error
import urllib.request


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--source", required=True, dest="url",
                   help="local TSV(.gz) from PARSE_SWISSPROT, or a URL serving the same TSV")
    p.add_argument("--out", required=True, help="Output JSON path")
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def fetch(url, out_path):
    def _download(ctx):
        req = urllib.request.Request(url, headers={"User-Agent": "domainsplit-pipeline"})
        with urllib.request.urlopen(req, context=ctx, timeout=600) as resp, open(out_path, "wb") as fh:
            while True:
                chunk = resp.read(1024 * 1024)
                if not chunk:
                    break
                fh.write(chunk)

    if url.startswith(("http://", "https://", "ftp://", "file://")):
        try:
            _download(ssl.create_default_context())
        except (urllib.error.URLError, ssl.SSLError) as exc:
            print(f"WARNING: SSL validation failed for {url} ({exc!r}); retrying unverified.",
                  file=sys.stderr, flush=True)
            _download(ssl._create_unverified_context())
    elif os.path.exists(url):
        shutil.copy(url, out_path)
    else:
        raise SystemExit(f"--source '{url}' is neither a URL nor a local file")


def open_maybe_gzip(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "rt")


def main():
    args = parse_args()

    raw = "swissprot.tsv"
    fetch(args.url, raw)

    accession_to_pfams = {}

    n_lines = 0
    with open_maybe_gzip(raw) as fh:
        for i, line in enumerate(fh):
            line = line.rstrip("\n")
            if i == 0 and line.lower().startswith("entry"):
                continue  # header
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                cols += [""] * (4 - len(cols))
            accession, pfam_field = cols[0], cols[3]
            if not accession:
                continue
            n_lines += 1

            pfams = sorted({p for p in pfam_field.replace(",", ";").split(";") if p})
            accession_to_pfams[accession] = pfams

    n_single = sum(1 for pfams in accession_to_pfams.values() if len(pfams) == 1)
    print(f"[swissprot_map] proteins={n_lines} single_domain={n_single}", flush=True)

    with open(args.out, "w") as fh:
        json.dump({"accession_to_pfams": accession_to_pfams}, fh)

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")


if __name__ == "__main__":
    main()
