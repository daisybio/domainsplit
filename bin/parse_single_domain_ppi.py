#!/usr/bin/env python3
"""Infer positive DDIs from HIPPIE PPIs between two single-domain proteins.

A PPI contributes a positive DDI only when *both* interactors are reviewed human
proteins annotated with exactly one Pfam domain; the DDI is then the pair of
those two single domains.  Identifiers in the HIPPIE columns may be UniProt
accessions or entry names (e.g. ``AL1A1_HUMAN``) -- both are resolved via the
SwissProt map.

Emits the normalized external-DDI TSV (``negative=0, source='single_domain_ppi'``);
insertion happens later, in ``insert_external_sources.py``.

HIPPIE ships in two layouts and the difference is silent, not loud: the classic
one is ``name_A, entrez_A, name_B, entrez_B, score, info`` while
``HIPPIE-current.txt`` is ``uniprot_accession_A, uniprot_name_A, entrez_A,
uniprot_accession_B, uniprot_name_B, entrez_B, score, info``. Reading the current
file at the classic offsets puts an entry name where the score is expected, every
row fails to parse as a float, and the step reports zero DDIs from a 170 MB input
without erroring. So the header decides the layout, and a file without one falls
back to the classic offsets.
"""

import argparse
import itertools
import json
import sys

from external_ddi_tsv import write_external_tsv


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--hippie", required=True)
    p.add_argument("--swissprot-map", required=True)
    p.add_argument("--min-score", type=float, required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


# Column offsets per layout: (identifier A candidates, identifier B candidates, score).
LEGACY_LAYOUT = ((0,), (2,), 4)
CURRENT_LAYOUT = ((0, 1), (3, 4), 6)


def detect_layout(first_line):
    """``(layout, is_header)`` for HIPPIE's two column layouts."""
    cols = [c.strip().lower() for c in first_line.rstrip("\n").split("\t")]
    if "score" in cols:
        a = tuple(i for i, c in enumerate(cols) if c.startswith("uniprot") and c.endswith("_a"))
        b = tuple(i for i, c in enumerate(cols) if c.startswith("uniprot") and c.endswith("_b"))
        layout = (a or CURRENT_LAYOUT[0], b or CURRENT_LAYOUT[1], cols.index("score"))
        return layout, True
    return LEGACY_LAYOUT, False


def main():
    args = parse_args()

    with open(args.swissprot_map) as fh:
        smap = json.load(fh)
    accession_to_pfams = smap["accession_to_pfams"]
    name_to_accession = smap["name_to_accession"]

    # single-domain proteins: accession -> its one Pfam
    single_domain = {
        acc: pfams[0]
        for acc, pfams in accession_to_pfams.items()
        if len(pfams) == 1
    }
    print(f"[single_domain_ppi] single-domain proteins: {len(single_domain)}", flush=True)

    def resolve_pfam(token):
        """Return the single Pfam of ``token`` (accession or name), else None."""
        acc = token if token in accession_to_pfams else name_to_accession.get(token)
        if acc is None:
            return None
        return single_domain.get(acc)

    def resolve_any(cols, indices):
        """First of ``indices`` whose identifier is a known single-domain protein."""
        for i in indices:
            if i < len(cols):
                pfam = resolve_pfam(cols[i])
                if pfam:
                    return pfam
        return None

    pairs = []
    n_rows = n_kept = n_unresolved = n_unscored = 0
    with open(args.hippie) as fh:
        first = fh.readline()
        (col_a, col_b, col_score), is_header = detect_layout(first)
        print(f"[single_domain_ppi] columns: A={col_a} B={col_b} score={col_score} "
              f"(header {'present' if is_header else 'absent'})", flush=True)
        lines = fh if is_header else itertools.chain([first], fh)
        for line in lines:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) <= col_score:
                continue
            n_rows += 1
            try:
                score = float(cols[col_score])
            except ValueError:
                n_unscored += 1
                continue
            if score < args.min_score:
                continue
            pfam_a = resolve_any(cols, col_a)
            pfam_b = resolve_any(cols, col_b)
            if pfam_a is None or pfam_b is None:
                n_unresolved += 1
                continue
            pairs.append((pfam_a, pfam_b))
            n_kept += 1
    if n_rows and n_unscored == n_rows:
        raise SystemExit(
            f"{args.hippie}: no row has a numeric score in column {col_score} -- "
            "the column layout was detected wrongly, so no DDI would be reported"
        )

    print(f"[single_domain_ppi] hippie_rows={n_rows} score>= {args.min_score}: "
          f"single_domain_pairs={n_kept} unresolved_or_multi={n_unresolved}", flush=True)

    rows = [(a, b, 0, "single_domain_ppi") for a, b in pairs]
    n = write_external_tsv(args.out, rows)
    print(f"[single_domain_ppi] wrote {n} rows to {args.out}", flush=True)

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")


if __name__ == "__main__":
    main()
