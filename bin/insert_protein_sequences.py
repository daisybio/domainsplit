#!/usr/bin/env python3
"""Fill in the sequence of every protein that carries a domain instance.

Reads a gzipped UniProt FASTA, filters it to the accessions named in the
protein-domain map, and upserts the sequence onto the ``protein`` row
INGEST_INSTANCES already created.

No embeddings. Protein-level per-residue vectors are retired: the models see the
*cut domain sequence* instead, and those embeddings are published as HDF5 files
(see ``export_domain_embeddings.py``) rather than pickled into SQLite. After this
change the master database carries no BLOBs at all.
"""

import argparse
import gzip
import shutil
import sqlite3
import sys

import pandas as pd
from Bio import SeqIO
from tqdm import tqdm


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-in", required=True)
    parser.add_argument("--uniprot-db", required=True)
    parser.add_argument("--protein-domain-map", required=True)
    parser.add_argument("--versions", required=True)
    parser.add_argument("--process-name", required=True)
    args = parser.parse_args()

    shutil.copy(args.db_in, "domainsplit.sqlite3")
    conn = sqlite3.connect("domainsplit.sqlite3")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    print("Inserting required uniprot records", flush=True)
    with (
        gzip.open(args.uniprot_db, "rt") as uniprot_text,
        gzip.open(args.protein_domain_map, "rt") as protein_domain_map_text,
    ):
        # Only `uniprot_id` is read: the map's `sequence` column holds domain
        # sequences, which INGEST_INSTANCES already wrote to domain_protein_map.
        pd_map = pd.read_csv(protein_domain_map_text, usecols=["uniprot_id"])
        required_uniprot_ids = set(pd_map["uniprot_id"].unique())

        uniprot_records = SeqIO.parse(uniprot_text, "fasta")
        uniprot_records = map(
            lambda record: (record.id.split("|")[1], str(record.seq)),
            uniprot_records,
        )
        uniprot_records = filter(
            lambda record: record[0] in required_uniprot_ids, uniprot_records
        )
        uniprot_records = tqdm(uniprot_records)

        # An upsert, not a plain insert: INGEST_INSTANCES has already created a
        # bare `protein` row for every parent protein of a domain instance, so a
        # plain INSERT would hit UNIQUE(uniprot_id) on the first record. This step
        # owns the sequence and fills it in.
        conn.executemany(
            """INSERT INTO protein (uniprot_id, sequence)
            VALUES (?, ?)
            ON CONFLICT(uniprot_id) DO UPDATE SET
                sequence = excluded.sequence;""",
            uniprot_records,
        )
    conn.commit()

    # Fail loudly on the one hole this step can leave. A protein carries a domain
    # instance only because FETCH_DOMAIN_META found it in the UniProt flat files;
    # if the FASTA carved out of those same files does not have it, the two inputs
    # were derived from different universes -- and the symptom otherwise is a NULL
    # sequence, no GO terms and nothing reporting either. Both now come from the
    # same `uniprot_dat_urls` list, so this cannot drift; it fails here rather than
    # being trusted not to.
    missing = conn.execute(
        "SELECT COUNT(*), MIN(uniprot_id) FROM protein WHERE sequence IS NULL"
    ).fetchone()
    total = conn.execute("SELECT COUNT(*) FROM protein").fetchone()[0]
    conn.close()
    print(f"[sequences] {total - missing[0]} of {total} proteins have a sequence", flush=True)
    if missing[0]:
        raise SystemExit(
            f"ERROR: {missing[0]} of {total} `protein` rows have no sequence "
            f"(e.g. {missing[1]}). Every parent protein of a domain instance came out "
            "of the UniProt flat files, so the FASTA parsed from those same files must "
            "cover it -- these two were built from different universes. Check that "
            "PARSE_SWISSPROT and FETCH_DOMAIN_META received the same "
            "--uniprot_dat_urls."
        )

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    pandas: {pd.__version__}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
