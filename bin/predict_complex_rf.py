#!/usr/bin/env python3

import argparse
import os
import sys
import shutil

import Bio
from Bio.PDB.PDBParser import PDBParser

import utils_struct
import pandas as pd

# Chain assignment: protein_a -> chain A, protein_b -> chain B.
RF2_CHAIN_A = "A"
RF2_CHAIN_B = "B"



def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--db_in", required=True)
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




def create_fasta(protein_id_a, protein_id_b, sequence_a, sequence_b, fastatmp_dir):
    fasta_header_a = f">{protein_id_a}\n"
    fasta_header_b = f">{protein_id_b}\n"
    fasta_filename = os.path.join(fastatmp_dir, f"{protein_id_a}_{protein_id_b}.fasta")
    with open(fasta_filename, "w") as fasta_file:
        fasta_file.write(fasta_header_a)
        fasta_file.write(sequence_a + "\n")
        fasta_file.write(fasta_header_b)
        fasta_file.write(sequence_b + "\n")

    return fasta_filename





def run_rf2_for_pair(protein_id_a: str, protein_id_b: str, sequence_a: str, sequence_b: str, outdir: str) -> None:
    """
    Invoke RosettaFold2 for the two sequences.

    TODO: replace stub with real RF2 CLI call, e.g.:
        rf_fold.py --fasta <fasta_file> --outdir <outdir>
    sequence_a must be passed as the first chain so it lands on chain A.
    """
    fasta_file = create_fasta(protein_id_a, protein_id_b, sequence_a, sequence_b, outdir)
    utils_struct.mock_predict_complex_rf(sequence_a, sequence_b, outdir)
    
    
    # pass   # ← implement RF2 invocation here


def _write_versions(versions_path: str, process_name: str) -> None:
    with open(versions_path, "w") as fh:
        fh.write(f'"{process_name}":\n')
        fh.write(f"    python: {sys.version.split()[0]}\n")
        fh.write(f"    biopython: {Bio.__version__}\n")
        fh.write(f"    RosettaFold2: TODO\n")



def subset_ppis(conn, ppis, limit: int):
    # Limit PPIs, such that for each DDI in the dataset, only the first N PPIs are kept.
    # This is useful for testing, to avoid running AF3 on all PPIs.
    if limit is None:
        return ppis

    print(f"[predict_complex_af] Limiting PPIs to {limit} per DDI", flush=True)

    query_ddi = "SELECT domain_id_a, domain_id_b FROM domain_domain_interaction"
    query_pd = "SELECT protein_id, domain_id FROM domain_protein_map"

    ddi = pd.read_sql(query_ddi, conn)
    pd_map = pd.read_sql(query_pd, conn)

    new_ppis = []

    for _, row in ddi.iterrows():
        domain_a = row["domain_id_a"]
        domain_b = row["domain_id_b"]

        proteins_a = pd_map[pd_map["domain_id"] == domain_a]["protein_id"].tolist()
        proteins_b = pd_map[pd_map["domain_id"] == domain_b]["protein_id"].tolist()

        ppis_subset = []
        for protein_a in proteins_a:
            for protein_b in proteins_b:
                ppis_subset.append((protein_a, protein_b))

        # Keep only the first N PPIs for this DDI
        ppis_subset = ppis_subset[:limit]
        # Add back the limited subset of PPIs for this DDI
        new_ppis.extend(ppis_subset)

    # Subset the original ppis list to only include the limited PPIs
    limited_ppis = [ppi for ppi in ppis if (ppi[0], ppi[3]) in new_ppis or (ppi[3], ppi[0]) in new_ppis]

    return limited_ppis


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    conn = utils_struct.connect_db(args.db_in)

    ppis = utils_struct.get_ppis(conn)
    limited_ppis = subset_ppis(conn, ppis, limit=5)  # Limit to 5 PPIs per DDI for testing
    ppis = limited_ppis
    conn.close()

    print(f"[predict_complex_rf] {len(ppis)} PPIs to process", flush=True)

    prediction_dir = os.path.join(args.outdir, "rf_predictions")
    os.makedirs(prediction_dir, exist_ok=True)

    for (protein_id_a, uniprot_id_a, seq_a, protein_id_b, uniprot_id_b, seq_b) in ppis:
        if not seq_a or not seq_b:
            print(f"  SKIP {uniprot_id_a}/{uniprot_id_b}: missing sequence", flush=True)
            continue

        # 1. Run RF2
        pair_dir = os.path.join(prediction_dir, f"{protein_id_a}_{protein_id_b}")
        os.makedirs(pair_dir, exist_ok=True)
        run_rf2_for_pair(protein_id_a, protein_id_b, seq_a, seq_b, pair_dir)

        # 2. Select best model
        pdb_path = select_rf_model(pair_dir)
        if pdb_path is None:
            print(f"  SKIP {protein_id_a}/{protein_id_b}: no RF2 model found", flush=True)
            continue
        # 3. Validate PDB
        # if not validate_pdb(pdb_path):
        #     print(f"  SKIP {protein_id_a}/{protein_id_b}: RF2 model failed validation", flush=True)
        #     continue

        # 4. Save PDB in outdir
        out_pdb_path = os.path.join(args.outdir, f"{protein_id_a}_{protein_id_b}.pdb")
        shutil.copy(pdb_path, out_pdb_path)

    # Remove the prediction_dir after processing all PPIs to save space
    shutil.rmtree(prediction_dir)

    _write_versions(args.versions, args.process_name)
    print("[predict_complex_rf] done", flush=True)



if __name__ == "__main__":
    main()