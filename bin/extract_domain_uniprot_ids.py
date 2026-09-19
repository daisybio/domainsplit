#!/usr/bin/env python3
"""Extract unique UniProt IDs from a domain FASTA and write a minimal FASTA
keyed by UniProt ID, for fetch_pdb_structures.py's alignment lookup.

Domain FASTA headers are '{family}_{uniprot_id}_{start}_{end}'
(fetch_domains.py). fetch_pdb_structures.py only reads record IDs off its
input FASTA (SeqIO.parse -> rec.id) -- it never looks at the sequence -- so
the placeholder sequence written here is never used for anything. One row
per distinct UniProt ID is all that's needed, not one per domain instance:
the same protein can appear behind dozens of domain headers, and there is
no reason to re-resolve/re-query the same UniProt ID that many times.
"""

import argparse
import gzip
import sys


def _open_fasta(path: str, mode: str):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-fasta", required=True, help="domain FASTA, keyed by {family}_{uniprot_id}_{start}_{end}")
    parser.add_argument("--output-fasta", required=True)
    args = parser.parse_args()

    from Bio import SeqIO

    uniprot_ids = set()
    skipped = 0
    with _open_fasta(args.input_fasta, "r") as fh:
        for rec in SeqIO.parse(fh, "fasta"):
            parts = rec.id.rsplit("_", 3)
            if len(parts) != 4:
                skipped += 1
                continue
            uniprot_ids.add(parts[1])

    if skipped:
        print(f"warn: {skipped} header(s) didn't parse as family_uniprot_start_end; skipped", file=sys.stderr)

    with open(args.output_fasta, "w") as out:
        for uid in sorted(uniprot_ids):
            out.write(f">{uid}\nX\n")

    print(f"{len(uniprot_ids)} unique UniProt IDs written to {args.output_fasta}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
