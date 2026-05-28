#!/usr/bin/env python3
"""Insert protein-protein interactions from STRING.

Reads the STRING PPI file (space-separated) and a gzipped UniProt ID mapping
to translate STRING protein IDs to UniProt IDs, then inserts PPI scores into
the protein_protein_interaction table.
"""

import argparse
import gzip
import shutil
import sqlite3
import sys
from typing import Dict

import pandas as pd
from tqdm import tqdm


def load_uniprot_id_mapping(mapping_path: str, key_name: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with gzip.open(mapping_path, "rt") as f:
        for line in f:
            uniprot_id, id_type, symbol_name = line.strip().split("\t")
            if id_type.strip() == key_name:
                mapping[symbol_name.strip()] = uniprot_id.strip()
    return mapping


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-in", required=True)
    parser.add_argument("--string-ppi", required=True)
    parser.add_argument("--uniprot-id-mapping", required=True)
    parser.add_argument("--versions", required=True)
    parser.add_argument("--process-name", required=True)
    args = parser.parse_args()

    shutil.copy(args.db_in, "domainsplit.sqlite3")
    conn = sqlite3.connect("domainsplit.sqlite3")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    print("Inserting PPI", flush=True)

    uniprot_to_pid = dict(
        conn.execute("SELECT uniprot_id, id FROM protein").fetchall()
    )
    print(f"Loaded {len(uniprot_to_pid)} protein ID mappings", flush=True)

    string_id_mapping = load_uniprot_id_mapping(args.uniprot_id_mapping, "STRING")
    ppi_df = pd.read_csv(args.string_ppi, sep=" ")

    insert_rows = []
    for _, row in tqdm(ppi_df.iterrows(), total=len(ppi_df)):
        uniprot_a = string_id_mapping.get(row["protein1"])
        uniprot_b = string_id_mapping.get(row["protein2"])
        if uniprot_a and uniprot_b:
            pid_a = uniprot_to_pid.get(uniprot_a)
            pid_b = uniprot_to_pid.get(uniprot_b)
            if pid_a is not None and pid_b is not None:
                insert_rows.append((pid_a, pid_b, row["combined_score"]))

    conn.executemany(
        "INSERT INTO protein_protein_interaction(protein_id_a, protein_id_b, score) "
        "VALUES (?, ?, ?)",
        insert_rows,
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
