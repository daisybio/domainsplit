#!/usr/bin/env python3
"""Insert domain-protein mapping with per-domain ESM embeddings.

Reads a gzipped protein-domain map CSV and looks up ESM3/ESMC per-domain
embeddings from an HDF5 file keyed as {pfam_id}_{uniprot_id}_{start}_{end}.
Inserts into the domain_protein_map table.
"""

import argparse
import gzip
import shutil
import sqlite3
import sys

import h5py
import numpy as np
import pandas as pd
from tqdm import tqdm


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-in", required=True)
    parser.add_argument("--protein-domain-map", required=True)
    parser.add_argument("--esm-domain-embeddings", required=True)
    parser.add_argument("--versions", required=True)
    parser.add_argument("--process-name", required=True)
    args = parser.parse_args()

    shutil.copy(args.db_in, "domainsplit.sqlite3")
    conn = sqlite3.connect("domainsplit.sqlite3")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    print("Inserting domain-protein mapping", flush=True)

    def iterate_domain_protein_alignments():
        with (
            gzip.open(args.protein_domain_map, "rt") as protein_domain_map_text,
            h5py.File(args.esm_domain_embeddings, mode="r") as domain_embeddings,
        ):
            pd_map = pd.read_csv(protein_domain_map_text)
            for row in pd_map.itertuples():
                key = f"{row.pfam_id}_{row.uniprot_id}_{row.start_pos}_{row.end_pos}"
                esm3_key = f"{key}/esm3"
                esmc_key = f"{key}/esmc"
                esm3_per_domain = esmc_per_domain = None
                if esm3_key in domain_embeddings:
                    esm3_per_domain = np.array(domain_embeddings[esm3_key]).dumps()
                if esmc_key in domain_embeddings:
                    esmc_per_domain = np.array(domain_embeddings[esmc_key]).dumps()
                yield (
                    row.sequence,
                    row.start_pos,
                    row.end_pos,
                    esm3_per_domain,
                    esmc_per_domain,
                    row.pfam_id,
                    row.uniprot_id,
                )

    conn.executemany(
        """
        INSERT OR IGNORE INTO domain_protein_map(
            domain_id, protein_id, domain_sequence,
            start_pos, end_pos,
            esm3_per_domain, esmc_per_domain
        )
        SELECT domain.id as domain_id, protein.id as protein_id,
            ? as domain_sequence, ? as start_pos, ? as end_pos,
            ? as esm3_per_domain, ? as esmc_per_domain
        FROM domain, protein
        WHERE
            domain.pfam_id = ? AND
            protein.uniprot_id = ?;
        """,
        tqdm(iterate_domain_protein_alignments()),
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
