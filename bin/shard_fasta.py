#!/usr/bin/env python3
"""Split a (gzipped) FASTA into N contiguous gzipped shards.

Used by SHARD_FASTA to fan out per-residue / pooled ESM embedding generation
across multiple GPU jobs. Contiguous (not round-robin) so each shard preserves
the input order — downstream length-bucketing still re-sorts internally, so
contiguous is simpler and faster to write.
"""

import argparse
import gzip
import sys
from pathlib import Path

from Bio import SeqIO


def _open_fasta(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-fasta", required=True, help="Path to input FASTA (.fasta or .fasta.gz)")
    parser.add_argument("--output-prefix", required=True, help="Output shard prefix (shards: <prefix>_<idx>.fasta.gz)")
    parser.add_argument("--num-shards", type=int, required=True)
    args = parser.parse_args()

    if args.num_shards < 1:
        sys.exit("--num-shards must be >= 1")

    with _open_fasta(args.input_fasta) as fh:
        records = list(SeqIO.parse(fh, "fasta"))

    total = len(records)
    if total == 0:
        sys.exit(f"No records in {args.input_fasta}")

    n = min(args.num_shards, total)
    chunk_size = (total + n - 1) // n

    for i in range(n):
        start = i * chunk_size
        end = min(start + chunk_size, total)
        if start >= end:
            break
        shard_path = Path(f"{args.output_prefix}_{i:03d}.fasta.gz")
        with gzip.open(shard_path, "wt") as out_fh:
            SeqIO.write(records[start:end], out_fh, "fasta")
        print(f"shard {i:03d}: {end - start} records -> {shard_path}", flush=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
