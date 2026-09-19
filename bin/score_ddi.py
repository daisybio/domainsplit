#!/usr/bin/env python3
"""
score_ddi.py
------------
Phase 2 of DDI scoring: apply ONE method's train-set empirical potential to
every DDI that method places in any split, and write two per-method files:

    <method>.scores.tsv    — ddi_id, method, split, instance_id_a,
                              instance_id_b, z_score
    <method>.confirmed.tsv — ddi_id, method, majority_confirmed, mean_confirmed

merge_scores.py applies these into the shared master once every method's
SCORE_DDI has run.

Scoring methodology (Sprinzak & Margalit 2001, PNAS):
    S_ab  = log10(C_ab / E_ab)                         if C_ab > 0 and E_ab >= 5
          = -0.5 (sc_mc) / -1.0 (sc_sc) / -0.5 (mc_mc)  otherwise
    E_ab  = T_DB * freq_a * freq_b
    z     = (score - mean_random) / std_random   over 1000 random trials
A domain pair is only scored if it has >= 5 interacting residue pairs.
"""
import argparse
import csv
import math
import random
import sys
from pathlib import Path

import h5py
import numpy as np
from Bio.PDB.PDBParser import PDBParser

from utils_struct import (
    bytes_to_tempfile, BACKBONE_ATOMS, CC_VDW, NO_SB,
    calculate_rsa_residue_level, extract_domain_residues,
    get_structure_bytes, get_structures_for_method,
)

TRIALS = 1000
MIN_INTERACTING_PAIRS = 5
C_CUTOFF = 0
E_CUTOFF = 5.0

T_DB: int = 0
DB_FREQ: dict = {}
C_AB_MATRIX: dict = {}
FALLBACK_SCORES = {'sc_sc': -1.0, 'sc_mc': -0.5, 'mc_mc': -0.5}
CHAIN_A, CHAIN_B = "A", "B"
ZSCORE_THRESHOLD = 2.3


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db_in", required=True)
    p.add_argument("--structures_h5", required=True)
    p.add_argument("--method", required=True)
    p.add_argument("--c_ab_matrix", required=True)
    p.add_argument("--db_freq", required=True)
    p.add_argument("--t_db", required=True)
    p.add_argument("--scores_out", required=True)
    p.add_argument("--confirmed_out", required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process_name", required=True)
    return p.parse_args()


def load_c_ab_matrix(path):
    global C_AB_MATRIX
    with open(path) as f:
        for row in csv.DictReader(f):
            C_AB_MATRIX[(row['aaA'], row['aaB'])] = float(row['count'])


def load_db_freq(path):
    global DB_FREQ
    with open(path) as f:
        for row in csv.DictReader(f):
            DB_FREQ[row['aa']] = float(row['frequency'])


def load_t_db(path):
    global T_DB
    T_DB = int(Path(path).read_text().strip())


def get_atom_interaction(atomA, atomB):
    elemA = (atomA.element or "").strip().upper()
    elemB = (atomB.element or "").strip().upper()
    pair = tuple(sorted([elemA, elemB]))
    dist = atomA - atomB
    qualifies = (
        (pair == ('C', 'C') and dist <= CC_VDW) or
        (pair == ('N', 'O') and dist <= NO_SB)
    )
    if not qualifies:
        return math.inf, None
    mc_a = atomA.get_id() in BACKBONE_ATOMS
    mc_b = atomB.get_id() in BACKBONE_ATOMS
    itype = 'mc_mc' if (mc_a and mc_b) else ('sc_mc' if (mc_a or mc_b) else 'sc_sc')
    return dist, itype


def check_for_interaction(resA, resB):
    best_dist, best_itype = math.inf, None
    for atomA in resA.get_atoms():
        for atomB in resB.get_atoms():
            dist, itype = get_atom_interaction(atomA, atomB)
            if dist < best_dist:
                best_dist, best_itype = dist, itype
    return best_itype


def check_rsa(domain, chain_id, threshold=0.1):
    rsa = calculate_rsa_residue_level(domain, chain_id=chain_id)
    return [rid for rid, v in rsa.items() if v is not None and v > threshold]


def find_interacting_residues(res_a, res_b, structure):
    n_a, n_b = len(res_a), len(res_b)
    matrix = [[""] * n_b for _ in range(n_a)]
    n_interacting = 0
    checked_a = check_rsa(structure, chain_id=CHAIN_A)
    checked_b = check_rsa(structure, chain_id=CHAIN_B)
    for i, rA in enumerate(res_a):
        if rA.get_id() not in checked_a:
            continue
        for j, rB in enumerate(res_b):
            if rB.get_id() not in checked_b:
                continue
            itype = check_for_interaction(rA, rB)
            if itype is not None:
                matrix[i][j] = itype
                n_interacting += 1
    return matrix, n_interacting


def get_expected_contacts(aaA, aaB):
    return T_DB * DB_FREQ.get(aaA, 0.0) * DB_FREQ.get(aaB, 0.0)


def compute_res_empirical_score(aaA, aaB, interaction_type):
    C_ab = C_AB_MATRIX.get((aaA, aaB), C_AB_MATRIX.get((aaB, aaA), 0.0))
    E_ab = get_expected_contacts(aaA, aaB)
    if C_ab > C_CUTOFF and E_ab >= E_CUTOFF:
        return math.log10(C_ab / E_ab)
    return FALLBACK_SCORES.get(interaction_type, -1.0)


def compute_empirical_potential(residues_a, residues_b, matrix):
    total = 0.0
    for i, rA in enumerate(residues_a):
        for j, rB in enumerate(residues_b):
            itype = matrix[i][j]
            if itype:
                aaA = rA if isinstance(rA, str) else rA.get_resname()
                aaB = rB if isinstance(rB, str) else rB.get_resname()
                total += compute_res_empirical_score(aaA, aaB, itype)
    return total


def sample_sequence(length, freq_dict):
    aas, weights = list(freq_dict.keys()), list(freq_dict.values())
    if sum(weights) == 0:
        weights = [1.0] * len(aas)
    return random.choices(aas, weights=weights, k=length)


def run_random_trials(n_a, n_b, matrix, n_trials):
    return [compute_empirical_potential(sample_sequence(n_a, DB_FREQ), sample_sequence(n_b, DB_FREQ), matrix)
            for _ in range(n_trials)]


def compute_z_score(real_score, random_scores):
    mean_r, std_r = float(np.mean(random_scores)), float(np.std(random_scores))
    return 0.0 if std_r == 0 else (real_score - mean_r) / std_r


def aggregate_scores(scores):
    aggregated = {}
    for ddi_id, score_list in scores.items():
        confirmed_votes = sum(c for _, c in score_list)
        mean_z = np.mean([z for z, _ in score_list])
        aggregated[ddi_id] = {
            'majority_confirmed': int(confirmed_votes > len(score_list) / 2),
            'mean_confirmed': int(mean_z >= ZSCORE_THRESHOLD),
        }
    return aggregated


def write_scores_tsv(path, method, rows):
    with open(path, "w") as f:
        f.write("ddi_id\tmethod\tsplit\tinstance_id_a\tinstance_id_b\tz_score\n")
        for ddi_id, split, ia, ib, z in rows:
            f.write(f"{ddi_id}\t{method}\t{split}\t{ia}\t{ib}\t{z}\n")


def write_confirmed_tsv(path, method, aggregated):
    with open(path, "w") as f:
        f.write("ddi_id\tmethod\tmajority_confirmed\tmean_confirmed\n")
        for ddi_id, s in aggregated.items():
            f.write(f"{ddi_id}\t{method}\t{s['majority_confirmed']}\t{s['mean_confirmed']}\n")


def write_versions(path, process_name):
    import Bio
    with open(path, "w") as f:
        f.write(f'"{process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    biopython: {Bio.__version__}\n")
        f.write(f"    numpy: {np.__version__}\n")
        f.write(f"    h5py: {h5py.__version__}\n")


def main():
    args = parse_args()
    print(f"[score_ddi] method={args.method}, zscore_threshold={ZSCORE_THRESHOLD}, trials={TRIALS}", flush=True)

    load_c_ab_matrix(args.c_ab_matrix)
    load_db_freq(args.db_freq)
    load_t_db(args.t_db)
    print(f"  Loaded T_DB={T_DB}, {len(C_AB_MATRIX)} C_ab entries, {len(DB_FREQ)} frequency entries", flush=True)

    method_structures = get_structures_for_method(args.db_in, args.method)
    print(f"[score_ddi] method={args.method}: {len(method_structures)} structures to score", flush=True)

    if not method_structures:
        print(f"[score_ddi] WARNING: no structures for method={args.method}, writing empty outputs", flush=True)
        write_scores_tsv(args.scores_out, args.method, [])
        write_confirmed_tsv(args.confirmed_out, args.method, {})
        write_versions(args.versions, args.process_name)
        return

    scores, score_rows = {}, []
    with h5py.File(args.structures_h5, "r") as h5file:
        for (ddi_id, instance_id_a, instance_id_b, split) in method_structures:
            pdb_gz = get_structure_bytes(h5file, instance_id_a, instance_id_b)
            if pdb_gz is None:
                print(f"WARNING: no structure bytes for ddi {ddi_id} ({instance_id_a}, {instance_id_b})", file=sys.stderr)
                continue

            path_to_tmp_pdb = bytes_to_tempfile(pdb_gz)
            structure = None
            try:
                structure = PDBParser(QUIET=True).get_structure(f"ddi_{ddi_id}", path_to_tmp_pdb)
            except Exception as e:
                print(f"WARNING: Failed to parse PDB for DDI {ddi_id}: {e}", file=sys.stderr)
            if structure is None:
                continue

            res_a = extract_domain_residues(structure, CHAIN_A)
            res_b = extract_domain_residues(structure, CHAIN_B)
            if not res_a or not res_b:
                print(f"  WARNING: no residues extracted for ddi {ddi_id}", file=sys.stderr)
                continue

            matrix, n_interacting = find_interacting_residues(res_a, res_b, structure)
            if n_interacting < MIN_INTERACTING_PAIRS:
                print(f"  Row {ddi_id}: only {n_interacting} interacting pairs "
                      f"(< {MIN_INTERACTING_PAIRS}) -- skipping", flush=True)
                continue

            real_score = compute_empirical_potential(res_a, res_b, matrix)
            random_scores = run_random_trials(len(res_a), len(res_b), matrix, TRIALS)
            z_score = compute_z_score(real_score, random_scores)
            confirmed = int(z_score >= ZSCORE_THRESHOLD)

            print(f"DDI {ddi_id}: n_interacting={n_interacting}, "
                  f"score={real_score:.3f}, z={z_score:.3f}, confirmed={confirmed}", flush=True)

            score_rows.append((ddi_id, split, instance_id_a, instance_id_b, z_score))
            scores.setdefault(ddi_id, []).append((z_score, confirmed))

    aggregated = aggregate_scores(scores)
    print(f"[score_ddi] method={args.method}: scored {len(score_rows)} instances, "
          f"{len(aggregated)} DDIs voted on", flush=True)

    write_scores_tsv(args.scores_out, args.method, score_rows)
    write_confirmed_tsv(args.confirmed_out, args.method, aggregated)
    write_versions(args.versions, args.process_name)
    print(f"[score_ddi] method={args.method}: done", flush=True)


if __name__ == "__main__":
    main()