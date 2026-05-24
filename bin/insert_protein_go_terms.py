#!/usr/bin/env python3
"""Insert protein GO term annotations from a UniProt GO terms TSV.

Reads a gzipped UniProt TSV with 'Entry' and 'Gene Ontology IDs' columns,
splits multi-valued GO IDs, and inserts into the protein_go_terms table by
joining on the protein table's uniprot_id.
"""

import argparse
import gzip
import itertools
import shutil
import sqlite3
import sys

import pandas as pd
from tqdm import tqdm


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-in", required=True)
    parser.add_argument("--uniprot-go-terms", required=True)
    parser.add_argument("--versions", required=True)
    parser.add_argument("--process-name", required=True)
    args = parser.parse_args()

    shutil.copy(args.db_in, "domainsplit.sqlite3")
    conn = sqlite3.connect("domainsplit.sqlite3")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    print("Inserting protein GO information", flush=True)
    with gzip.open(args.uniprot_go_terms, "rt") as go_terms_text:
        go_terms_df = pd.read_csv(go_terms_text, sep="\t")
        insert_iterator = (
            (
                (go_term, row["Entry"])
                for go_term in row["Gene Ontology IDs"].split("; ")
            )
            for _, row in go_terms_df.iterrows()
            if isinstance(row["Gene Ontology IDs"], str)
        )
        insert_iterator = itertools.chain.from_iterable(insert_iterator)
        insert_iterator = tqdm(insert_iterator)

        conn.executemany(
            """
            INSERT INTO protein_go_terms(protein_id, go_accession)
            SELECT protein.id as protein_id, ? AS go_accession
            FROM protein
            WHERE protein.uniprot_id = ?;
            """,
            insert_iterator,
        )
    conn.commit()
    conn.close()

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    pandas: {pd.__version__}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
