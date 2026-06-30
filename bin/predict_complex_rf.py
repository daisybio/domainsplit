#! /usr/bin/env python3

import argparse
import os
import shutil
import sqlite3
import sys

import Bio
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO

import utils_struct

# Chain assignment: protein_a -> chain A, protein_b -> chain B.
RF2_CHAIN_A = "A"
RF2_CHAIN_B = "B"



def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--db_in", required=True)
    p.add_argument("--db_out", required=True)
    p.add_argument("--outdir", required=True, help="Directory for PPI PDB files")
    p.add_argument("--versions", required=True)
    p.add_argument("--fastatmp", required=True, help="Directory to store tmp fasta files")
    p.add_argument("--process_name", required=True)
    return p.parse_args()


def select_rf_model(model_dir: str):
    """Return path to the best RF2 model (model_0.pdb), or None."""
    candidate = os.path.join(model_dir, "model_0.pdb")
    return candidate if os.path.exists(candidate) else None


def validate_pdb(pdb_path: str) -> bool:
    try:
        PDBParser(QUIET=True).get_structure("model", pdb_path)
        return True
    except Exception as exc:
        print(f"  WARNING: PDB validation failed for {pdb_path}: {exc}", flush=True)
        return False


def create_fasta():
    pass





def run_rf2_for_pair(uniprot_id_a: str, uniprot_id_b: str, sequence_a: str, sequence_b: str, outdir: str) -> None:
    """
    Invoke RosettaFold2 for the two sequences.

    TODO: replace stub with real RF2 CLI call, e.g.:
        rf_fold.py --fasta <fasta_file> --outdir <outdir>
    sequence_a must be passed as the first chain so it lands on chain A.
    """
    
    pass   # ← implement RF2 invocation here


def _write_versions(versions_path: str, process_name: str) -> None:
    with open(versions_path, "w") as fh:
        fh.write(f'"{process_name}":\n')
        fh.write(f"    python: {sys.version.split()[0]}\n")
        fh.write(f"    biopython: {Bio.__version__}\n")
        fh.write(f"    RosettaFold2: TODO\n")




def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    shutil.copy(args.db_in, args.db_out)
    conn = utils_struct.connect_db(args.db_out)
    utils_struct.create_tables(conn)

    ppis = utils_struct.get_ppis(conn)
    print(f"[predict_complex_rf] {len(ppis)} PPIs to process", flush=True)

    prediction_dir = os.path.join(args.outdir, "rf_predictions")
    os.makedirs(prediction_dir, exist_ok=True)

    for (protein_id_a, uniprot_id_a, seq_a, protein_id_b, uniprot_id_b, seq_b) in ppis:
        if not seq_a or not seq_b:
            print(f"  SKIP {uniprot_id_a}/{uniprot_id_b}: missing sequence", flush=True)
            continue

        # 1. Run RF2
        pair_dir = os.path.join(prediction_dir, f"{uniprot_id_a}_{uniprot_id_b}")
        os.makedirs(pair_dir, exist_ok=True)
        run_rf2_for_pair(uniprot_id_a, uniprot_id_b, seq_a, seq_b, pair_dir)

        pdb_path = select_rf_model(pair_dir)
        if pdb_path is None:
            print(f"  SKIP {uniprot_id_a}/{uniprot_id_b}: no RF2 model found", flush=True)
            continue
        if not validate_pdb(pdb_path):
            print(f"  SKIP {uniprot_id_a}/{uniprot_id_b}: RF2 model failed validation", flush=True)
            continue

        # 2. Record chain assignment
        utils_struct.store_chain_map(conn, protein_id_a, protein_id_b, RF2_CHAIN_A, RF2_CHAIN_B, 'RF2')
        conn.commit()

    conn.close()
    _write_versions(args.versions, args.process_name)
    print("[predict_complex_rf] done", flush=True)



if __name__ == "__main__":
    main()