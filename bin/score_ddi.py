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

PERFORMANCE NOTE (vs. original):
  - Contact detection (which residue pairs interact, and at what mc/sc
    type) is vectorized with NumPy/cdist instead of a Python double-loop
    over every atom pair. This is normally >90% of the per-structure
    runtime.
  - The 1000 random trials no longer re-scan the full n_a x n_b matrix in
    Python each time. A global (20,20,3) score lookup table is built once,
    and all 1000 trials are sampled and scored via vectorized NumPy
    indexing in one shot.
  - Semantics are preserved exactly: for each trial a single random AA is
    drawn per residue *position* (not per contact), matching the original
    sample_sequence() behaviour where a residue's identity is shared
    across all of its contacts within one trial.
"""
import argparse
import csv
import math
import sys
from pathlib import Path

import h5py
import numpy as np
from scipy.spatial.distance import cdist
from Bio.PDB.PDBParser import PDBParser

from utils_struct import (
    bytes_to_tempfile, BACKBONE_ATOMS, CC_VDW, NO_SB, AA_3, AA_SET,
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
# Relative solvent accessibility of the unbound protein must be >= 10% to
# count a residue as surface (excludes buried side-chains from scoring).
RSA_THRESHOLD = 0.1

# Fixed ordering for itype indices used throughout the vectorized code.
ITYPES = ['mc_mc', 'sc_mc', 'sc_sc']
ITYPE_INDEX = {t: i for i, t in enumerate(ITYPES)}

AA_INDEX = {aa: i for i, aa in enumerate(AA_3)}

# Populated once in main() after C_AB_MATRIX / DB_FREQ / T_DB are loaded.
global SCORE_LOOKUP
SCORE_LOOKUP = np.array([])  # shape (20, 20, 3), see build_score_lookup()


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


def get_expected_contacts(aaA, aaB):
    return T_DB * DB_FREQ.get(aaA, 0.0) * DB_FREQ.get(aaB, 0.0)


def compute_res_empirical_score(aaA, aaB, interaction_type):
    C_ab = C_AB_MATRIX.get((aaA, aaB), C_AB_MATRIX.get((aaB, aaA), 0.0))
    E_ab = get_expected_contacts(aaA, aaB)
    if C_ab > C_CUTOFF and E_ab >= E_CUTOFF:
        return math.log10(C_ab / E_ab)
    return FALLBACK_SCORES.get(interaction_type, -1.0)


def build_score_lookup():
    """
    Build the (20, 20, 3) empirical-potential lookup table once, globally.
    table[i, j, k] = compute_res_empirical_score(AA_3[i], AA_3[j], ITYPES[k])
    Only 20*20*3 = 1200 combinations total -- computed once instead of being
    re-derived per residue pair per trial.
    """
    table = np.zeros((len(AA_3), len(AA_3), len(ITYPES)))
    for i, aaA in enumerate(AA_3):
        for j, aaB in enumerate(AA_3):
            for k, itype in enumerate(ITYPES):
                table[i, j, k] = compute_res_empirical_score(aaA, aaB, itype)
    return table


def check_rsa(domain, chain_id, threshold=RSA_THRESHOLD):
    """
    Residue ids with relative solvent accessibility >= threshold (10% by
    default, per Sprinzak & Margalit 2001: buried side-chains excluded).
    """
    rsa = calculate_rsa_residue_level(domain, chain_id=chain_id)
    return {rid for rid, v in rsa.items() if v is not None and v >= threshold}


def _residue_atom_arrays(residues, checked_ids):
    """
    Flatten C/N/O atoms of the given residues (restricted to those in
    checked_ids, i.e. RSA-passing) into parallel NumPy arrays:
        coords     (n_atoms, 3) float
        elems      (n_atoms,)   'C'/'N'/'O'
        is_backbone(n_atoms,)   bool
        res_idx    (n_atoms,)   int index into `residues`
    Only C/N/O atoms are kept since those are the only elements the
    original contact definition (C-C vdW, N-O salt bridge) ever uses.
    """
    coords, elems, is_backbone, res_idx = [], [], [], []
    for i, r in enumerate(residues):
        if r.get_id() not in checked_ids:
            continue
        for atom in r.get_atoms():
            el = (atom.element or "").strip().upper()
            if el not in ("C", "N", "O"):
                continue
            coords.append(atom.coord)
            elems.append(el)
            is_backbone.append(atom.get_id() in BACKBONE_ATOMS)
            res_idx.append(i)
    return (
        np.asarray(coords, dtype=float).reshape(-1, 3),
        np.asarray(elems),
        np.asarray(is_backbone, dtype=bool),
        np.asarray(res_idx, dtype=int),
    )


def find_contacts_vectorized(res_a, res_b, structure):
    """
    Vectorized replacement for find_interacting_residues/check_for_interaction.

    Returns (ra_arr, rb_arr, itype_idx_arr, n_interacting):
      - ra_arr, rb_arr: residue indices (into res_a / res_b) for each
        contacting residue pair
      - itype_idx_arr: index into ITYPES for that pair's contact type,
        taken from the closest qualifying atom pair (matches the original
        "best_dist" semantics)
      - n_interacting: number of contacting residue pairs
    """
    checked_a = check_rsa(structure, chain_id=CHAIN_A)
    checked_b = check_rsa(structure, chain_id=CHAIN_B)

    coords_a, elems_a, mc_a, ridx_a = _residue_atom_arrays(res_a, checked_a)
    coords_b, elems_b, mc_b, ridx_b = _residue_atom_arrays(res_b, checked_b)

    n_b = len(res_b)
    empty = (np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=int), 0)
    if len(coords_a) == 0 or len(coords_b) == 0:
        return empty

    dist = cdist(coords_a, coords_b)

    cc_mask = (elems_a[:, None] == "C") & (elems_b[None, :] == "C") & (dist <= CC_VDW)
    no_mask = (
        ((elems_a[:, None] == "N") & (elems_b[None, :] == "O")) |
        ((elems_a[:, None] == "O") & (elems_b[None, :] == "N"))
    ) & (dist <= NO_SB)
    hit = cc_mask | no_mask
    if not hit.any():
        return empty

    ai, bi = np.where(hit)
    d = dist[ai, bi]
    ra = ridx_a[ai]
    rb = ridx_b[bi]
    mc_pair_a = mc_a[ai]
    mc_pair_b = mc_b[bi]

    # itype per atom-pair hit: mc_mc if both backbone, sc_sc if neither, else sc_mc
    itype_idx = np.where(
        mc_pair_a & mc_pair_b, ITYPE_INDEX['mc_mc'],
        np.where(~mc_pair_a & ~mc_pair_b, ITYPE_INDEX['sc_sc'], ITYPE_INDEX['sc_mc'])
    )

    # Reduce to one entry per (residue_a, residue_b) pair: the closest
    # qualifying atom-pair, matching the original best_dist behaviour.
    key = ra.astype(np.int64) * n_b + rb.astype(np.int64)
    order = np.argsort(d)
    key_sorted = key[order]
    itype_sorted = itype_idx[order]
    _, first_idx = np.unique(key_sorted, return_index=True)

    best_key = key_sorted[first_idx]
    best_itype = itype_sorted[first_idx]
    best_ra = (best_key // n_b).astype(int)
    best_rb = (best_key % n_b).astype(int)

    return best_ra, best_rb, best_itype, len(best_key)


def compute_real_score_vectorized(res_a, res_b, ra_arr, rb_arr, itype_arr):
    """Empirical potential for the actual (non-random) residue identities."""
    if len(ra_arr) == 0:
        return 0.0
    aa_a_idx = np.array([AA_INDEX.get(res_a[i].get_resname(), -1) for i in ra_arr])
    aa_b_idx = np.array([AA_INDEX.get(res_b[j].get_resname(), -1) for j in rb_arr])
    valid = (aa_a_idx >= 0) & (aa_b_idx >= 0)
    if not valid.any():
        return 0.0
    return float(SCORE_LOOKUP[aa_a_idx[valid], aa_b_idx[valid], itype_arr[valid]].sum())


def run_random_trials_vectorized(n_a, n_b, ra_arr, rb_arr, itype_arr, n_trials):
    """
    Vectorized replacement for run_random_trials/compute_empirical_potential.

    A random AA identity is drawn per residue *position* (length n_a and
    n_b respectively) for each of n_trials trials -- exactly matching the
    original sample_sequence() semantics, where a residue's random
    identity is shared across all of its contacts within one trial.
    Contact-level scores are then gathered via NumPy fancy indexing
    instead of a per-trial Python loop over the interaction matrix.
    """
    if len(ra_arr) == 0:
        return np.zeros(n_trials)

    aas = list(DB_FREQ.keys())
    weights = np.array([DB_FREQ[a] for a in aas], dtype=float)
    weights = weights / weights.sum() if weights.sum() > 0 else np.full(len(aas), 1.0 / len(aas))
    aa_pool_idx = np.array([AA_INDEX.get(a, -1) for a in aas])
    valid_pool = aa_pool_idx >= 0
    aa_pool_idx = aa_pool_idx[valid_pool]
    weights = weights[valid_pool]
    weights = weights / weights.sum()

    sample_a_full = np.random.choice(aa_pool_idx, size=(n_trials, n_a), p=weights)
    sample_b_full = np.random.choice(aa_pool_idx, size=(n_trials, n_b), p=weights)

    aa_a_per_contact = sample_a_full[:, ra_arr]  # (n_trials, n_contacts)
    aa_b_per_contact = sample_b_full[:, rb_arr]  # (n_trials, n_contacts)

    scores = SCORE_LOOKUP[aa_a_per_contact, aa_b_per_contact, itype_arr[None, :]]
    return scores.sum(axis=1)  # (n_trials,)


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
    import scipy
    with open(path, "w") as f:
        f.write(f'"{process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    biopython: {Bio.__version__}\n")
        f.write(f"    numpy: {np.__version__}\n")
        f.write(f"    scipy: {scipy.__version__}\n")
        f.write(f"    h5py: {h5py.__version__}\n")


def main():
    global SCORE_LOOKUP
    args = parse_args()
    print(f"[score_ddi] method={args.method}, zscore_threshold={ZSCORE_THRESHOLD}, trials={TRIALS}", flush=True)

    load_c_ab_matrix(args.c_ab_matrix)
    load_db_freq(args.db_freq)
    load_t_db(args.t_db)
    print(f"  Loaded T_DB={T_DB}, {len(C_AB_MATRIX)} C_ab entries, {len(DB_FREQ)} frequency entries", flush=True)

    SCORE_LOOKUP = build_score_lookup()
    print("  Built global (20,20,3) score lookup table", flush=True)

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

            ra_arr, rb_arr, itype_arr, n_interacting = find_contacts_vectorized(res_a, res_b, structure)
            if n_interacting < MIN_INTERACTING_PAIRS:
                print(f"  Row {ddi_id}: only {n_interacting} interacting pairs "
                      f"(< {MIN_INTERACTING_PAIRS}) -- skipping", flush=True)
                continue

            real_score = compute_real_score_vectorized(res_a, res_b, ra_arr, rb_arr, itype_arr)
            random_scores = run_random_trials_vectorized(
                len(res_a), len(res_b), ra_arr, rb_arr, itype_arr, TRIALS
            )
            z_score = compute_z_score(real_score, random_scores)
            confirmed = int(z_score >= ZSCORE_THRESHOLD)

            # print(f"DDI {ddi_id}: n_interacting={n_interacting}, "
            #       f"score={real_score:.3f}, z={z_score:.3f}, confirmed={confirmed}", flush=True)

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