#!/usr/bin/env python3

import argparse
import os
import sys
# import shutil

import Bio
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO

import utils_struct
import pandas as pd


# Chain assignment: protein_a -> chain A, protein_b -> chain B.
AF3_CHAIN_A = "A"
AF3_CHAIN_B = "B"

path_to_af_script_monomer = "/nfs/data/alphafold3/scripts/run_monomer_data_pipelines.sh"
path_to_af_script_complex = "/nfs/data/alphafold3/scripts/run_pairwise_inferences.sh"

path_to_input_data = "/nfs/data/alphafold3/input/c.thomas/"
path_to_output_data = "/nfs/data/alphafold3/output/c.thomas/"
cpus = "8"
dataset_name = "teschting_smoke_100_dataset"
out_name = "teschting_smoke_100_batch"



def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--db_in", required=True)
    p.add_argument("--outdir", required=True, help="Directory for PPI PDB files")
    p.add_argument("--versions", required=True)
    p.add_argument("--process_name", required=True)
    return p.parse_args()



def select_af_models(model_dir: str, file_ending: str):
    """Return path to the best AF3 model (*.cif), or None."""
    # Each ppi has a directory with the following structure:
    # <model_dir>/
    # ├── <protein_id_a>_<protein_id_b>/
    #     ├── <process_id>/
    #         ├── <protein_id_a>_<protein_id_b>/
    #           ├── <protein_id_a>_<protein_id_b>_model.cif

    pair_results = {}

    for dir in os.listdir(model_dir):
        pair_dir = os.path.join(model_dir, dir)
        pair_results[dir] = ""
        if not os.path.isdir(pair_dir):
            continue

        for process_id in os.listdir(pair_dir):
            process_dir = os.path.join(pair_dir, process_id)
            if not os.path.isdir(process_dir):
                continue

            model_subdir = os.path.join(process_dir, dir)
            if not os.path.isdir(model_subdir):
                continue

            for file in os.listdir(model_subdir):
                if file.endswith(file_ending):
                    pair_results[dir] = os.path.join(model_subdir, file)
                    break

    return pair_results


def convert_cif_to_pdb(cif_path: str, pdb_path: str) -> None:
    try:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("model", cif_path)
        io_obj = PDBIO()
        io_obj.set_structure(structure)
        io_obj.save(pdb_path)
    except Exception as e:
        print(f"Error converting {cif_path} to {pdb_path}: {e}", flush=True)


def prep_af_input(ppi_data, dir) -> None:
    """
    Write input tables for AF3 to disk, in the format expected by AF3.
    1. tsv file with protein_1  protein_2, called pairs.tsv
    2. tsv file with protein_id	sequence, called proteins.tsv
    """
    protein_sequences = {}
    pairs = []
    for protein_id_a, uniprot_id_a, seq_a, protein_id_b, uniprot_id_b, seq_b in ppi_data:
        protein_sequences[protein_id_a] = seq_a
        protein_sequences[protein_id_b] = seq_b
        pairs.append((protein_id_a, protein_id_b))

    pairs_path = os.path.join(dir, "pairs.tsv")
    proteins_path = os.path.join(dir, "proteins.tsv")
    
    pd.DataFrame(pairs, columns=["protein_1", "protein_2"]).to_csv(pairs_path, sep="\t", index=False)
    pd.DataFrame(protein_sequences.items(), columns=["protein_id", "sequence"]).to_csv(proteins_path, sep="\t", index=False)



def _write_versions(versions_path: str, process_name: str) -> None:
    with open(versions_path, "w") as fh:
        fh.write(f'"{process_name}":\n')
        fh.write(f"    python: {sys.version.split()[0]}\n")
        fh.write(f"    biopython: {Bio.__version__}\n")


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
    limited_ppis = subset_ppis(conn, ppis, limit=2)  # Limit to 2 PPIs per DDI for testing
    ppis = limited_ppis
    conn.close()

    print(f"[predict_complex_af] {len(ppis)} PPIs to process", flush=True)

    # prediction_dir = os.path.join(args.outdir, "af_predictions")
    # os.makedirs(prediction_dir, exist_ok=True)

    dataset_dir = os.path.join(path_to_input_data, dataset_name)
    os.makedirs(dataset_dir, exist_ok=True)

    dataset_dir_out = os.path.join(path_to_output_data, out_name)
    os.makedirs(dataset_dir_out, exist_ok=True)

    # Create input files for AF3
    prep_af_input(ppis, dataset_dir)

    pair_inference_dir = os.path.join(dataset_dir_out, "af3_outputs", "pair_inferences")

    all_models = select_af_models(pair_inference_dir, "cif")

    # For each found model, copy it to the output directory with the name protein_id_a_protein_id_b.pdb (after converting from .cif to .pdb)
    for name, cif_path in all_models.items():
        if cif_path:
            protein_id_a, protein_id_b = name.split("__")
            pdb_path = os.path.join(args.outdir, f"{protein_id_a}_{protein_id_b}.pdb")
            convert_cif_to_pdb(cif_path, pdb_path)
        else:
            print(f"  SKIP {name}: no AF3 model found", flush=True)
    
    _write_versions(args.versions, args.process_name)
    print("[predict_complex_af] done", flush=True)


if __name__ == "__main__":
    main()