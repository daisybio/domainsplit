#!/usr/bin/env python3
"""Parse Negatome into the normalized external-DDI TSV.

Negatome ``combined_pfam.txt`` lists whitespace-separated Pfam pairs that do not
interact.  Emitted as ``negative=1, source='negatome'``; nothing is written to
the database here -- external sources are inserted last, by
``insert_external_sources.py``, so that 3did and the sampled negatives keep
first claim on a pair.

Each line must start with two Pfam accessions (a ``.version`` suffix is
stripped).  Anything else is counted and sampled, not passed on: the download
falls back to an unverified TLS context and the default host is an old MIPS
server, so an HTML error page served with HTTP 200 is a realistic input, and
read token-by-token it used to become a handful of ``<!DOCTYPE``/``html>``
"domains" and a green run with no Negatome DDIs.  A file with no valid pair at
all is fatal for the same reason.
"""

import argparse
import re
import sys

from external_ddi_tsv import write_external_tsv


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--negatome", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


PFAM = re.compile(r"(PF\d{5})(?:\.\d+)?")


def iter_negatome_pairs(path, rejected=None):
    """Yield ``(pfam_a, pfam_b)``; append every other non-blank line to ``rejected``."""
    with open(path) as f:
        for line in f:
            tokens = line.split()
            if not tokens:
                continue
            if len(tokens) >= 2:
                a, b = PFAM.fullmatch(tokens[0]), PFAM.fullmatch(tokens[1])
                if a and b:
                    yield a.group(1), b.group(1)
                    continue
            if rejected is not None:
                rejected.append(line.rstrip("\n"))


def main():
    args = parse_args()

    rejected = []
    pairs = list(iter_negatome_pairs(args.negatome, rejected))
    print(f"[negatome] read {len(pairs)} pairs", flush=True)
    if rejected:
        sample = "; ".join(repr(line[:80]) for line in rejected[:5])
        print(f"[negatome] WARNING: skipped {len(rejected)} lines that are not two Pfam "
              f"accessions, e.g. {sample}", flush=True)
    if not pairs:
        raise SystemExit(
            f"ERROR: {args.negatome} holds no Pfam pair ({len(rejected)} non-blank lines, "
            "none starting with two PFxxxxx accessions). The download most likely returned "
            "an error page -- check --url_negatome."
        )

    rows = [(a, b, 1, "negatome") for a, b in pairs]
    n = write_external_tsv(args.out, rows)
    print(f"[negatome] wrote {n} rows to {args.out}", flush=True)

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")


if __name__ == "__main__":
    main()
