#!/usr/bin/env python3
"""Build a reviewed-human UniProt -> Pfam map for single-domain detection.

Downloads one UniProt stream (TSV, fields accession,id,gene_names,xref_pfam for
``reviewed:true AND organism_id:9606``) and emits ``swissprot_pfam_map.json``:

    {
      "accession_to_pfams": {accession: [Pfam, ...]},
      "name_to_accession":  {entry_name_or_gene: accession}
    }

``name_to_accession`` lets the single-domain step resolve HIPPIE identifiers that
are entry names (e.g. ``AL1A1_HUMAN``) or gene names; accessions resolve directly
against ``accession_to_pfams``.  Gene names that map to more than one accession are
dropped as ambiguous; unique entry names always win.
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
    p.add_argument("--url", required=True, help="UniProt stream URL or local TSV(.gz) file")
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
        raise SystemExit(f"url_uniprot_swissprot_pfam '{url}' is neither a URL nor a local file")


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
    gene_to_accs = {}      # gene token -> set of accessions (for ambiguity check)
    entry_name_to_acc = {}

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
            accession, entry_name, gene_names, pfam_field = cols[0], cols[1], cols[2], cols[3]
            if not accession:
                continue
            n_lines += 1

            pfams = sorted({p for p in pfam_field.replace(",", ";").split(";") if p})
            accession_to_pfams[accession] = pfams

            if entry_name:
                entry_name_to_acc[entry_name] = accession
            for token in gene_names.split():
                gene_to_accs.setdefault(token, set()).add(accession)

    # entry names are unique and authoritative; add unambiguous gene names that
    # do not collide with an entry name
    name_to_accession = dict(entry_name_to_acc)
    for token, accs in gene_to_accs.items():
        if token in name_to_accession:
            continue
        if len(accs) == 1:
            name_to_accession[token] = next(iter(accs))

    n_single = sum(1 for pfams in accession_to_pfams.values() if len(pfams) == 1)
    print(f"[swissprot_map] proteins={n_lines} single_domain={n_single} "
          f"names={len(name_to_accession)}", flush=True)

    with open(args.out, "w") as fh:
        json.dump(
            {"accession_to_pfams": accession_to_pfams,
             "name_to_accession": name_to_accession},
            fh,
        )

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")


if __name__ == "__main__":
    main()
