#!/usr/bin/env python3
"""
build_scoring_matrix.py
------------------------
Phase 1 of DDI scoring: derive ONE splitting method's empirical potential
from that method's own TRAIN-split positive-DDI structures (real AF3
predictions only -- mocked placeholders excluded).

Deliberately per-method: score_ddi.py scores every DDI under this method
(train/val/test alike) against the matrix built here. Built from anything
but this method's own train split, a test-time DDI would be scored against
a background that had already seen it.

Produces, per method:
    <method>.c_ab_matrix.csv — 20x20 contact counts (aaA, aaB, count)
    <method>.db_freq.csv     — train-set surface AA frequencies (aa, frequency)
    <method>.t_db.txt        — total contact count T_DB
"""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import Bio
import h5py
import numpy as np
from Bio.PDB.PDBParser import PDBParser

from utils_struct import (
    AA_3, bytes_to_tempfile, extract_domain_residues,
    get_positive_train_structures, get_structure_bytes, residues_contact,
)

CHAIN_A, CHAIN_B = "A", "B"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db_in", required=True)
    p.add_argument("--structures_h5", required=True)
    p.add_argument("--method", required=True)
    p.add_argument("--c_ab_matrix", required=True)
    p.add_argument("--db_freq", required=True)
    p.add_argument("--t_db", required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process_name", required=True)
    return p.parse_args()


def build_matrix(args):
    C_ab = defaultdict(lambda: defaultdict(float))
    surface_counts = defaultdict(int)
    n_instances = n_skipped = n_no_contact = n_missing_bytes = 0

    train_structures = get_positive_train_structures(args.db_in, args.method)
    print(f"[build_scoring_matrix] method={args.method}: {len(train_structures)} "
          f"train-split positive-DDI structures", flush=True)

    with h5py.File(args.structures_h5, "r") as h5file:
        for (ddi_id, instance_id_a, instance_id_b) in train_structures:
            n_instances += 1
            pdb_gz = get_structure_bytes(h5file, instance_id_a, instance_id_b)
            if pdb_gz is None:
                print(f"WARNING: no structure bytes for ddi {ddi_id} "
                      f"({instance_id_a}, {instance_id_b})", file=sys.stderr)
                n_missing_bytes += 1
                continue

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

            for r in res_a:
                surface_counts[r.get_resname()] += 1
            for r in res_b:
                surface_counts[r.get_resname()] += 1

            found_any = False
            for rA in res_a:
                aaA = rA.get_resname()
                for rB in res_b:
                    aaB = rB.get_resname()
                    if residues_contact(rA, rB):
                        C_ab[aaA][aaB] += 1
                        C_ab[aaB][aaA] += 1
                        found_any = True
            if not found_any:
                n_no_contact += 1
            if n_instances % 500 == 0:
                print(f"  Processed {n_instances} DDIs...", flush=True)

    print(f"Done. Instances: {n_instances}, skipped: {n_skipped}, "
          f"missing structure bytes: {n_missing_bytes}, no contacts: {n_no_contact}", flush=True)
    return C_ab, surface_counts


def write_c_ab_matrix(C_ab, path):
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["aaA", "aaB", "count"])
        for aaA in AA_3:
            for aaB in AA_3:
                writer.writerow([aaA, aaB, C_ab[aaA][aaB]])
    print(f"Written C_ab matrix to {path}", flush=True)


def write_db_freq(surface_counts, path):
    total = sum(surface_counts.values())
    if total == 0:
        print("WARNING: no surface residues accumulated — DB_FREQ will be zeros", file=sys.stderr)
        total = 1
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["aa", "frequency"])
        for aa in AA_3:
            writer.writerow([aa, surface_counts.get(aa, 0) / total])
    print(f"Written DB_FREQ to {path}", flush=True)


def write_t_db(C_ab, path):
    t_db_diag = int(sum(C_ab[aaA][aaA] for aaA in AA_3))
    t_db_offdiag = int(sum(
        C_ab[aaA][aaB] for i, aaA in enumerate(AA_3) for j, aaB in enumerate(AA_3) if j > i
    ))
    t_db = t_db_diag + 2 * t_db_offdiag
    Path(path).write_text(str(t_db) + "\n")
    print(f"T_DB = {t_db}, written to {path}", flush=True)
    return t_db


def write_versions(path, process_name):
    with open(path, "w") as f:
        f.write(f'"{process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    biopython: {Bio.__version__}\n")
        f.write(f"    numpy: {np.__version__}\n")
        f.write(f"    h5py: {h5py.__version__}\n")


def main():
    args = parse_args()
    print(f"[build_scoring_matrix] method={args.method}: building train-set scoring matrix", flush=True)
    C_ab, surface_counts = build_matrix(args)
    write_c_ab_matrix(C_ab, args.c_ab_matrix)
    write_db_freq(surface_counts, args.db_freq)
    t_db = write_t_db(C_ab, args.t_db)
    write_versions(args.versions, args.process_name)
    print(f"[build_scoring_matrix] method={args.method}: done, T_DB={t_db}", flush=True)


if __name__ == "__main__":
    main()