#! /usr/bin/env python3

import argparse
import os
import shutil
import sys

import Bio
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO

import utils_struct


# Chain assignment: protein_a -> chain A, protein_b -> chain B.
AF3_CHAIN_A = "A"
AF3_CHAIN_B = "B"



def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--db_in", required=True)
    p.add_argument("--db_out", required=True)
    p.add_argument("--outdir", required=True, help="Directory for PPI PDB files")
    p.add_argument("--versions", required=True)
    p.add_argument("--process_name", required=True)
    return p.parse_args()



def select_af_model(model_dir: str):
    """Return path to the best AF3 model (*_model_0.cif), or None."""
    for fname in os.listdir(model_dir):
        if fname.endswith("_model_0.cif"):
            return os.path.join(model_dir, fname)
    return None


def convert_cif_to_pdb(cif_path: str, pdb_path: str) -> None:
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("model", cif_path)
    io_obj = PDBIO()
    io_obj.set_structure(structure)
    io_obj.save(pdb_path)


def run_af3_for_pair(uniprot_id_a: str, uniprot_id_b: str, sequence_a: str, sequence_b: str, outdir: str) -> None:
    """
    Invoke AlphaFold 3 for the two sequences.

    TODO: replace stub with real AF3 CLI call, e.g.:
        af3 --seq1 <sequence_a> --seq2 <sequence_b> --outdir <outdir>
    sequence_a must be passed as the first chain so it lands on chain A.
    """
    utils_struct.mock_predict_complex(sequence_a, sequence_b, outdir, f"{uniprot_id_a}_{uniprot_id_b}")
    # pass   # ← implement AF3 invocation here



def _write_versions(versions_path: str, process_name: str) -> None:
    with open(versions_path, "w") as fh:
        fh.write(f'"{process_name}":\n')
        fh.write(f"    python: {sys.version.split()[0]}\n")
        fh.write(f"    biopython: {Bio.__version__}\n")
        fh.write(f"    AF3: TODO\n")




def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    shutil.copy(args.db_in, args.db_out)
    conn = utils_struct.connect_db(args.db_out)
    utils_struct.create_tables(conn)

    ppis = utils_struct.get_ppis(conn)
    print(f"[predict_complex_af] {len(ppis)} PPIs to process", flush=True)

    prediction_dir = os.path.join(args.outdir, "af_predictions")
    os.makedirs(prediction_dir, exist_ok=True)

    for (protein_id_a, uniprot_id_a, seq_a, protein_id_b, uniprot_id_b, seq_b) in ppis:
        if not seq_a or not seq_b:
            print(f"  SKIP {uniprot_id_a}/{uniprot_id_b}: missing sequence", flush=True)
            continue

        # 1. Run Alphafold3
        pair_dir = os.path.join(prediction_dir, f"{uniprot_id_a}_{uniprot_id_b}")
        os.makedirs(pair_dir, exist_ok=True)
        run_af3_for_pair(uniprot_id_a, uniprot_id_b, seq_a, seq_b, pair_dir)

        cif_path = select_af_model(pair_dir)
        if cif_path is None:
            print(f"  SKIP {uniprot_id_a}/{uniprot_id_b}: no AF3 model found", flush=True)
            continue

        # Convert cif to pdb
        pdb_path = os.path.join(pair_dir, f"{uniprot_id_a}_{uniprot_id_b}.pdb")
        # convert_cif_to_pdb(cif_path, pdb_path)

        # Record which uniprot_id id is associated with which chain
        # sequence_a → chain A, sequence_b → chain B
        utils_struct.store_chain_map(conn, protein_id_a, protein_id_b, AF3_CHAIN_A, AF3_CHAIN_B, 'AF3')

        conn.commit()

    conn.close()
    _write_versions(args.versions, args.process_name)
    print("[predict_complex_af] done", flush=True)


if __name__ == "__main__":
    main()