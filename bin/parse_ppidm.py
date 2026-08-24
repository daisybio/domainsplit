#!/usr/bin/env python3
"""Parse PPIDM predictions into the normalized external-DDI TSV.

Input ``predicted_ddi_ppi.tsv`` columns: ``domain_1  domain_2  class`` where each
domain token looks like ``10114/PF00069`` (the Pfam accession follows the slash).

Every kept DDI is emitted **twice**: once under the umbrella source ``PPIDM`` and
once under ``PPIDM_<Class>``.  Both land on the same database row, whose merged
source list reads ``"PPIDM,PPIDM_Gold"`` -- so a consumer can select all of PPIDM
or one confidence class without either query needing to know about the other.
Classes are emitted Gold -> Silver -> Bronze so the merged list is ordered by
confidence.
"""

import argparse
import sys

from external_ddi_tsv import write_external_tsv

# highest confidence first
CLASS_ORDER = ["Gold", "Silver", "Bronze"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--ppidm", required=True)
    p.add_argument("--classes", required=True,
                   help="comma-separated classes to include, e.g. 'Bronze,Silver,Gold'")
    p.add_argument("--out", required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def extract_pfam(token):
    """``10114/PF00069`` or ``PF00069.3`` -> ``PF00069`` (or None)."""
    pf = token.split("/")[-1].split(".")[0].strip()
    return pf if pf.startswith("PF") else None


def main():
    args = parse_args()

    allowed = {c.strip().capitalize() for c in args.classes.split(",") if c.strip()}
    classes = [c for c in CLASS_ORDER if c in allowed]
    print(f"[ppidm] including classes: {classes}", flush=True)

    pairs_by_class = {c: [] for c in classes}
    n_rows = n_bad = 0
    with open(args.ppidm) as fh:
        for i, line in enumerate(fh):
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                continue
            if i == 0 and cols[2].strip().lower() == "class":
                continue  # header
            cls = cols[2].strip().capitalize()
            if cls not in pairs_by_class:
                continue
            pfam_a = extract_pfam(cols[0])
            pfam_b = extract_pfam(cols[1])
            if pfam_a is None or pfam_b is None:
                n_bad += 1
                continue
            pairs_by_class[cls].append((pfam_a, pfam_b))
            n_rows += 1

    print(f"[ppidm] parsed {n_rows} pairs ({n_bad} unparseable)", flush=True)

    rows = []
    for cls in classes:  # Gold first
        for pfam_a, pfam_b in pairs_by_class[cls]:
            rows.append((pfam_a, pfam_b, 0, "PPIDM"))
            rows.append((pfam_a, pfam_b, 0, f"PPIDM_{cls}"))

    n = write_external_tsv(args.out, rows)
    print(f"[ppidm] wrote {n} rows to {args.out}", flush=True)

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")


if __name__ == "__main__":
    main()
