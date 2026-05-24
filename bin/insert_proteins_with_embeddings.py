#!/usr/bin/env python3
"""Insert proteins with their per-residue embeddings.

Reads a gzipped UniProt FASTA, filters to proteins referenced in the
protein-domain map, loads ProtT5 and ESM (esm3/esmc) per-residue embeddings
from HDF5 files, and inserts into the protein table.
"""

import argparse
import gzip
import shutil
import sqlite3
import sys

import h5py
import numpy as np
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-in", required=True)
    parser.add_argument("--uniprot-db", required=True)
    parser.add_argument("--protein-domain-map", required=True)
    parser.add_argument("--prott5-embeddings", required=True)
    parser.add_argument("--esm-protein-embeddings", required=True)
    parser.add_argument("--versions", required=True)
    parser.add_argument("--process-name", required=True)
    args = parser.parse_args()

    shutil.copy(args.db_in, "domainsplit.sqlite3")
    conn = sqlite3.connect("domainsplit.sqlite3")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    print("Inserting required uniprot records and embeddings", flush=True)
    with (
        gzip.open(args.uniprot_db, "rt") as uniprot_text,
        gzip.open(args.protein_domain_map, "rt") as protein_domain_map_text,
        h5py.File(args.esm_protein_embeddings, mode="r") as esm_protein_embeddings,
        h5py.File(args.prott5_embeddings, mode="r") as prott5_embeddings_file,
    ):
        pd_map = pd.read_csv(protein_domain_map_text)
        required_uniprot_ids = set(pd_map["uniprot_id"].unique())

        uniprot_records = SeqIO.parse(uniprot_text, "fasta")
        uniprot_records = map(
            lambda record: (record.id.split("|")[1], str(record.seq)),
            uniprot_records,
        )
        uniprot_records = filter(
            lambda record: record[0] in required_uniprot_ids, uniprot_records
        )

        def get_prott5_embedding(seq_id):
            if seq_id in prott5_embeddings_file:
                return np.array(prott5_embeddings_file[seq_id]).dumps()
            return None

        uniprot_records = map(
            lambda record: record + (get_prott5_embedding(record[0]),),
            uniprot_records,
        )

        def get_esm_embeddings(record):
            seq_id = record[0]
            esm3_key = f"{seq_id}/esm3"
            esmc_key = f"{seq_id}/esmc"
            esm3_embedding = (
                np.array(esm_protein_embeddings[esm3_key]).dumps()
                if esm3_key in esm_protein_embeddings
                else None
            )
            esmc_embedding = (
                np.array(esm_protein_embeddings[esmc_key]).dumps()
                if esmc_key in esm_protein_embeddings
                else None
            )
            return (esm3_embedding, esmc_embedding)

        uniprot_records = map(
            lambda record: record + get_esm_embeddings(record), uniprot_records
        )
        uniprot_records = tqdm(uniprot_records)

        conn.executemany(
            """INSERT INTO protein (
                uniprot_id, sequence,
                prott5_per_residue, esm3_per_residue, esmc_per_residue
            )
            VALUES (?, ?, ?, ?, ?);""",
            uniprot_records,
        )
    conn.commit()
    conn.close()

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    pandas: {pd.__version__}\n")
        f.write(f"    numpy: {np.__version__}\n")
        f.write(f"    h5py: {h5py.__version__}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
