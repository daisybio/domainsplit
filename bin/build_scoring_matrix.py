#!/usr/bin/env python3
"""
build_scoring_matrix.py
------------------------
Phase 1 of DDI scoring: derive the database-wide empirical potential
components from the full 3did flat file and all downloaded PDB structures.

This runs ONCE and produces three artefact files consumed by score_ddi.py:
    c_ab_matrix.csv  — 20×20 contact counts C_ab (rows: aaA, aaB, count)
    db_freq.csv      — database-wide surface AA frequencies (rows: aa, frequency)
    t_db.txt         — total contact count T_DB as a single integer

Methodology (from Sprinzak & Margalit 2001, PNAS):
    - For every #=3D instance in the 3did flat file:
        - Load the PDB file
        - Extract domain A residues (chain_a, res_start_a, res_end_a)
        - Extract domain B residues (chain_b, res_start_b, res_end_b)
        - For each residue pair (i, j):
            - Apply contact definition:
                H-bond:     N-O distance <= 3.5 Å
                Salt bridge: N-O distance <= 5.5 Å
                vdW:        C-C distance <= 5.0 Å
            - If contact found: accumulate C_ab[aaA][aaB] += 1
            - Record surface residue for DB_FREQ
    - T_DB = sum of all C_ab values
    - DB_FREQ[aa] = count of aa on domain surfaces / total surface residues
"""

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path
import Bio
import numpy as np

from Bio.PDB.PDBParser import PDBParser
from utils_struct import bytes_to_tempfile, AA_3, extract_domain_residues, residues_contact, get_domain_structures


CHAIN_A = "A"
CHAIN_B = "B"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--db_in", required=True,
                   help="Path to the SQLite database after Domain Slicing")
    p.add_argument("--c_ab_matrix", required=True,
                   help="Output CSV for C_ab matrix (aaA, aaB, count)")
    p.add_argument("--db_freq", required=True,
                   help="Output CSV for database surface AA frequencies")
    p.add_argument("--t_db", required=True,
                   help="Output TXT for total contact count T_DB")
    p.add_argument("--versions", required=True,
                   help="Output versions.yml")
    p.add_argument("--process_name", required=True,
                   help="Nextflow process name for versions.yml")
    return p.parse_args()






def build_matrix(args):
    """Build C_ab contact matrix and calculate surface residue frequencies"""

    # C_ab[aaA][aaB] — raw contact counts by amino acid type pair
    C_ab = defaultdict(lambda: defaultdict(float))

    # surface_counts[aa] — how many times this AA appeared on a domain surface
    surface_counts = defaultdict(int)

    n_instances  = 0
    n_skipped    = 0
    n_no_contact = 0


    for (ds_id, ddi_id, pdb_gz) in get_domain_structures(args.db_in):
        n_instances += 1
        
        path_to_tmp_pdb = bytes_to_tempfile(pdb_gz)
        structure = None
        try:
            structure = PDBParser(QUIET=True).get_structure(f"ddi_{ddi_id}", path_to_tmp_pdb)
        except Exception as e:
            print(f"WARNING: Failed to parse PDB for DDI {ddi_id}: {e}", file=sys.stderr)
        
        if structure is None:
            n_skipped += 1
            continue

        res_a = extract_domain_residues(structure, CHAIN_A)
        res_b = extract_domain_residues(structure, CHAIN_B)

        if not res_a or not res_b:
            n_skipped += 1
            continue

        # Accumulate surface residue frequencies
        for r in res_a:
            surface_counts[r.get_resname()] += 1
        for r in res_b:
            surface_counts[r.get_resname()] += 1

        # Count contacts by AA type pair
        found_any = False
        for rA in res_a:
            aaA = rA.get_resname()
            for rB in res_b:
                aaB = rB.get_resname()
                if residues_contact(rA, rB):
                    C_ab[aaA][aaB] += 1
                    C_ab[aaB][aaA] += 1  # symmetric
                    found_any = True

        if not found_any:
            n_no_contact += 1

        if n_instances % 500 == 0:
            print(f"  Processed {n_instances} DDIs...", flush=True)

    print(f"Done. Instances: {n_instances}, skipped: {n_skipped}, "
          f"no contacts found: {n_no_contact}", flush=True)

    return C_ab, surface_counts


def write_c_ab_matrix(C_ab, path):
    """Write C_ab as a CSV with columns: aaA, aaB, count."""
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["aaA", "aaB", "count"])
        for aaA in AA_3:
            for aaB in AA_3:
                count = C_ab[aaA][aaB]
                writer.writerow([aaA, aaB, count])
    print(f"Written C_ab matrix to {path}", flush=True)


def write_db_freq(surface_counts, path):
    """
    Write database-wide surface AA frequencies as CSV: aa, frequency.
    frequency = count / total_surface_residues  (sums to 1.0).
    """
    total = sum(surface_counts.values())
    if total == 0:
        print("WARNING: no surface residues accumulated — DB_FREQ will be zeros",
              file=sys.stderr)
        total = 1  # avoid division by zero

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["aa", "frequency"])
        for aa in AA_3:
            freq = surface_counts.get(aa, 0) / total
            writer.writerow([aa, freq])
    print(f"Written DB_FREQ to {path}", flush=True)


def write_t_db(C_ab, path):
    """Write T_DB = sum of all C_ab values´"""
    t_db_diag = int(sum(
        C_ab[aaA][aaA]
        for aaA in AA_3
    ))
    t_db_offdiag = int(sum(
        C_ab[aaA][aaB]        for i, aaA in enumerate(AA_3)
        for j, aaB in enumerate(AA_3) if j > i
    ))
    t_db = t_db_diag + 2 * t_db_offdiag  # off-diagonal counts are symmetric, so multiply by 2
    
    Path(path).write_text(str(t_db) + "\n")
    print(f"T_DB = {t_db}, written to {path}", flush=True)
    return t_db


def write_versions(path, process_name):
    
    python_version = sys.version.split()[0]
    with open(path, "w") as f:
        f.write(f'"{process_name}":\n')
        f.write(f"    python: {python_version}\n")
        f.write(f"    biopython: {Bio.__version__}\n")
        f.write(f"    numpy: {np.__version__}\n")


def main():
    args = parse_args()

    C_ab, surface_counts = build_matrix(args)

    write_c_ab_matrix(C_ab,         args.c_ab_matrix)
    write_db_freq(surface_counts,   args.db_freq)
    t_db = write_t_db(C_ab,         args.t_db)
    write_versions(args.versions,   args.process_name)

    print(f"Scoring matrix built successfully. T_DB = {t_db}", flush=True)


if __name__ == "__main__":
    main()