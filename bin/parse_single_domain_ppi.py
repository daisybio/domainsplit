#!/usr/bin/env python3
"""Infer positive DDIs from HIPPIE PPIs between two single-domain proteins.

A PPI contributes a positive DDI only when *both* interactors are reviewed human
proteins annotated with exactly one Pfam domain; the DDI is then the pair of
those two single domains.  Interactors are read from HIPPIE's **UniProt accession**
columns and nothing else: the accession is the identifier every other table in
this pipeline is keyed by, while an entry name (``AL1A1_HUMAN``) has to be
resolved through a second map, and gene names are ambiguous often enough that the
old resolver dropped them.

Emits the normalized external-DDI TSV (``negative=0, source='single_domain_ppi'``);
insertion happens later, in ``insert_external_sources.py``.

``HIPPIE-current.txt`` is ``uniprot_accession_A, uniprot_name_A, entrez_A,
uniprot_accession_B, uniprot_name_B, entrez_B, score, info``, and the accession
columns are located by name from the header. The older layout
(``name_A, entrez_A, name_B, entrez_B, score, info``) carries no accession at all
and is **rejected**, loudly: reading it at these offsets would put an entry name
where the score belongs, every row would fail to parse as a float, and the step
would report zero DDIs from a 170 MB input without erroring. A hard failure
naming the missing column is the only outcome that cannot be mistaken for "this
source contributed nothing".
"""

import argparse
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


#: Where the accession columns sit in HIPPIE-current.txt, used only if a header is
#: present but does not name them explicitly.
DEFAULT_ACCESSION_COLUMNS = (0, 3, 6)


def detect_layout(first_line):
    """``(col_a, col_b, col_score)`` from HIPPIE's header. Raises if it has none.

    Located by name rather than by offset: HIPPIE has changed its column set
    before, and a positional read of the wrong layout fails silently -- an entry
    name parses as a bad float, every row is skipped, and the step reports zero
    DDIs from a 170 MB file.
    """
    cols = [c.strip().lower() for c in first_line.rstrip("\n").split("\t")]
    if "score" not in cols:
        raise SystemExit(
            "ERROR: the HIPPIE file has no header row naming a 'score' column, so it "
            "is the older layout that carries no UniProt accessions. Supply "
            "HIPPIE-current.txt, whose columns are uniprot_accession_A, "
            "uniprot_name_A, entrez_A, uniprot_accession_B, uniprot_name_B, "
            "entrez_B, score, info."
        )
    a = [i for i, c in enumerate(cols) if c.startswith("uniprot_accession") and c.endswith("_a")]
    b = [i for i, c in enumerate(cols) if c.startswith("uniprot_accession") and c.endswith("_b")]
    if not a or not b:
        d_a, d_b, _ = DEFAULT_ACCESSION_COLUMNS
        a, b = a or [d_a], b or [d_b]
    return a[0], b[0], cols.index("score")


def main():
    args = parse_args()

    with open(args.swissprot_map) as fh:
        smap = json.load(fh)
    accession_to_pfams = smap["accession_to_pfams"]

    # single-domain proteins: accession -> its one Pfam
    single_domain = {
        acc: pfams[0]
        for acc, pfams in accession_to_pfams.items()
        if len(pfams) == 1
    }
    print(f"[single_domain_ppi] single-domain proteins: {len(single_domain)}", flush=True)

    def resolve_pfam(cols, index):
        """The single Pfam of the accession in column ``index``, else None."""
        if index >= len(cols):
            return None
        return single_domain.get(cols[index].strip())

    pairs = []
    n_rows = n_kept = n_unresolved = n_unscored = 0
    with open(args.hippie) as fh:
        first = fh.readline()
        col_a, col_b, col_score = detect_layout(first)
        print(
            f"[single_domain_ppi] accession columns: A={col_a} B={col_b} score={col_score}",
            flush=True,
        )
        # detect_layout only returns on a header row, so `first` is consumed.
        for line in fh:
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
            pfam_a = resolve_pfam(cols, col_a)
            pfam_b = resolve_pfam(cols, col_b)
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
