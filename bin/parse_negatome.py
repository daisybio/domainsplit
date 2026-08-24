#!/usr/bin/env python3
"""Parse Negatome into the normalized external-DDI TSV.

Negatome ``combined_pfam.txt`` lists whitespace-separated Pfam pairs that do not
interact.  Emitted as ``negative=1, source='negatome'``; nothing is written to
the database here -- external sources are inserted last, by
``insert_external_sources.py``, so that 3did and the sampled negatives keep
first claim on a pair.
"""

import argparse
import sys

from external_ddi_tsv import write_external_tsv


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--negatome", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def iter_negatome_pairs(path):
    with open(path) as f:
        for line in f:
            tokens = line.split()
            if len(tokens) < 2:
                continue
            yield tokens[0], tokens[1]


def main():
    args = parse_args()

    pairs = list(iter_negatome_pairs(args.negatome))
    print(f"[negatome] read {len(pairs)} pairs", flush=True)

    rows = [(a, b, 1, "negatome") for a, b in pairs]
    n = write_external_tsv(args.out, rows)
    print(f"[negatome] wrote {n} rows to {args.out}", flush=True)

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")


if __name__ == "__main__":
    main()
