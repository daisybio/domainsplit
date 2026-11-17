#!/usr/bin/env python3
"""Stream a BioGRID or IntAct MITAB release into a compact PPI TSV.

Output (gzip TSV): ``acc_a acc_b taxid method source_db score pmids``, one row
per kept input row. Accessions are suffix-stripped but not yet mapped to the
primary accession; that happens in build_ppidm_sources.py against UniProt.

The input is the release .zip (the relevant member is streamed, never extracted:
IntAct's is 11.6 GB) or a plain/gzipped MITAB text.
"""

import argparse
import fnmatch
import gzip
import io
import re
import sys
import zipfile
from collections import Counter

MI = re.compile(r"MI:(\d+)")
TAX = re.compile(r"taxid:(-?\d+)")
PUBMED = re.compile(r"pubmed:(\d+)")
MISCORE = re.compile(r"intact-miscore:([0-9.eE+-]+)")
BIOGRID_PHYSICAL = {"0915", "0407", "0403", "0914"}
MEMBER = {"biogrid": "*.mitab.txt", "intact": "intact.txt"}


def log(msg):
    print(f"[extract_ppi_mitab] {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True)
    p.add_argument("--format", required=True, choices=["biogrid", "intact"])
    p.add_argument("--out", required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def open_text(path, fmt):
    if zipfile.is_zipfile(path):
        zf = zipfile.ZipFile(path)
        members = [n for n in zf.namelist() if fnmatch.fnmatch(n.split("/")[-1], MEMBER[fmt])]
        if len(members) != 1:
            raise SystemExit(f"ERROR: {path}: expected exactly one member matching "
                             f"{MEMBER[fmt]!r}, found {members}")
        log(f"streaming {members[0]} from {path}")
        return io.TextIOWrapper(zf.open(members[0]), encoding="utf-8", errors="replace")
    with open(path, "rb") as fh:
        gz = fh.read(2) == b"\x1f\x8b"
    return (gzip.open if gz else open)(path, "rt", encoding="utf-8", errors="replace")


def strip_suffix(acc):
    """``P12345-2`` / ``P12345-PRO_0000012345`` -> ``P12345``."""
    return acc.split("-", 1)[0]


def mi(field):
    m = MI.search(field)
    return m.group(1) if m else ""


def taxid(field):
    m = TAX.search(field)
    return int(m.group(1)) if m else -1


def biogrid_side(field):
    """The single swiss-prot accession of a BioGRID alt-id field; None if absent, '' if ambiguous."""
    accs = {strip_suffix(e[len("uniprot/swiss-prot:"):])
            for e in field.split("|") if e.startswith("uniprot/swiss-prot:")}
    if not accs:
        return None
    return accs.pop() if len(accs) == 1 else ""


def main():
    args = parse_args()
    drops = Counter()
    n_rows = n_kept = 0
    with open_text(args.input, args.format) as fh, \
            gzip.open(args.out, "wt", compresslevel=6) as out:
        out.write("acc_a\tacc_b\ttaxid\tmethod\tsource_db\tscore\tpmids\n")
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 15:
                drops["short_row"] += 1
                continue
            n_rows += 1
            if n_rows % 1_000_000 == 0:
                log(f"{n_rows} rows read, {n_kept} kept")

            if args.format == "biogrid":
                if mi(c[11]) not in BIOGRID_PHYSICAL:
                    drops["non_physical_type"] += 1
                    continue
                a, b = biogrid_side(c[2]), biogrid_side(c[3])
                if a is None or b is None:
                    drops["no_swissprot_id"] += 1
                    continue
                if a == "" or b == "":
                    drops["several_swissprot_ids"] += 1
                    continue
                source_db, score = "biogrid", "-"
            else:
                if not (c[0].startswith("uniprotkb:") and c[1].startswith("uniprotkb:")):
                    drops["non_uniprot_interactor"] += 1
                    continue
                if len(c) > 35 and c[35].strip().lower() == "true":
                    drops["negative"] += 1
                    continue
                a = strip_suffix(c[0][len("uniprotkb:"):])
                b = strip_suffix(c[1][len("uniprotkb:"):])
                source_db = mi(c[12])
                scores = [float(s) for s in MISCORE.findall(c[14])]
                score = repr(max(scores)) if scores else "-"

            ta, tb = taxid(c[9]), taxid(c[10])
            if ta <= 0 or tb <= 0:
                drops["taxid_missing_or_negative"] += 1
                continue
            if ta != tb:
                drops["cross_taxon"] += 1
                continue
            if a == b:
                drops["self_pair"] += 1
                continue
            if a > b:
                a, b = b, a
            pmids = ";".join(sorted(set(PUBMED.findall(c[8]))))
            out.write(f"{a}\t{b}\t{ta}\t{mi(c[6])}\t{source_db}\t{score}\t{pmids}\n")
            n_kept += 1

    log(f"{args.format}: rows={n_rows} kept={n_kept} dropped={dict(sorted(drops.items()))}")
    if n_kept == 0:
        raise SystemExit(f"ERROR: no row of {args.input} survived; wrong --format or file?")

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")


if __name__ == "__main__":
    main()
