#!/usr/bin/env python3
"""Insert protein-protein interactions from STRING.

Reads the STRING PPI file (space-separated) and a gzipped UniProt ID mapping to
translate STRING protein IDs to UniProt IDs, then inserts PPI scores into the
protein_protein_interaction table.

The mapping now comes from PARSE_SWISSPROT's ``DR   STRING;`` lines rather than a
per-organism ``<ORG>_<taxid>_idmapping.dat.gz`` download; the format is the same
``uniprot \t id_type \t symbol`` triple, so nothing here had to change for it.
Nothing here is single-species either: STRING ids carry their own taxon prefix
(``9606.ENSP…``) and are globally unique, so a multi-species mapping and a
multi-species links file join correctly with no taxon bookkeeping at all.

Every dropped edge is counted. It used to be silent, which meant a mapping that
covered none of the links file -- the wrong organism, or a taxon whose entries the
``--taxon-ids`` filter excluded -- produced zero PPI rows and a green run.
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
    print(f"Loaded {len(string_id_mapping)} STRING -> UniProt mappings", flush=True)
    ppi_df = pd.read_csv(args.string_ppi, sep=" ")

    insert_rows = []
    n_unmapped = n_not_in_db = 0
    for _, row in tqdm(ppi_df.iterrows(), total=len(ppi_df)):
        uniprot_a = string_id_mapping.get(row["protein1"])
        uniprot_b = string_id_mapping.get(row["protein2"])
        if not uniprot_a or not uniprot_b:
            n_unmapped += 1
            continue
        pid_a = uniprot_to_pid.get(uniprot_a)
        pid_b = uniprot_to_pid.get(uniprot_b)
        if pid_a is None or pid_b is None:
            n_not_in_db += 1
            continue
        insert_rows.append((pid_a, pid_b, row["combined_score"]))

    # The two drops mean different things and only one of them is ordinary.
    # `not_in_db` is expected and usually the large number: most STRING proteins
    # carry no domain this run cares about. `unmapped` is the one to watch -- at
    # scale it means the links file and the mapping disagree about which species
    # they cover.
    print(
        f"STRING edges: {len(ppi_df)} read, {len(insert_rows)} inserted, "
        f"{n_unmapped} with an unmapped STRING id, {n_not_in_db} whose proteins "
        "are not in this database",
        flush=True,
    )
    if len(ppi_df) and not insert_rows:
        raise SystemExit(
            f"ERROR: none of {len(ppi_df)} STRING edges could be inserted "
            f"({n_unmapped} unmapped, {n_not_in_db} absent from the database). "
            "The STRING links file and the UniProt mapping most likely cover "
            "different species -- check --url_string against --swissprot_taxon_ids."
        )

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
