#!/usr/bin/env python3
"""
score_ddi.py
------------
Phase 2 of DDI scoring: apply the precomputed 3did empirical potential to
all predicted complex structures and update ddi_complex_coords in the DB.

Only rows where z_score IS NULL are processed -- PDB-sourced rows (already
scored by GET_PDB directly from the 3did flat file) are left untouched.

Inputs (produced once by build_scoring_matrix.py):
    --matrix    c_ab_matrix.csv  (20x20 contact counts, database-wide)
    --dbfreq    db_freq.csv      (database-wide surface AA frequencies)
    --t-db      t_db.txt         (total contact count T_DB)

Scoring methodology (Sprinzak & Margalit 2001, PNAS):
    S_ab  = log10(C_ab / E_ab)                         if C_ab > 0 and E_ab >= 5
          = -0.5 (sc_mc) / -1.0 (sc_sc) / -0.5 (mc_mc)  otherwise
    E_ab  = T_DB * freq_a * freq_b
    score = sum of S_ab over all interacting residue pairs in the complex
    z     = (score - mean_random) / std_random   over 1000 random trials

Contact definition (used both to decide "is this residue pair interacting"
and to classify the interaction as sc_sc / sc_mc / mc_mc):
    H-bond:      N-O distance <= 3.5 A
    Salt bridge: N-O distance <= 5.5 A
    vdW contact: C-C distance <= 5.0 A
A domain pair is only scored at all if it has >= 5 interacting residue pairs
(the 3did minimum-contact requirement).
"""


import argparse
import csv
import math
import random
import shutil
import sqlite3
import sys
from pathlib import Path

import numpy as np
from Bio.PDB.PDBParser import PDBParser
from utils_struct import bytes_to_tempfile, BACKBONE_ATOMS, CC_VDW, NO_SB, calculate_rsa_residue_level, extract_domain_residues, update_score, get_domain_structures

# Constants
TRIALS = 1000
MIN_INTERACTING_PAIRS = 5

# Cutoffs applied at SCORE LOOKUP time (per residue-type-pair), not when
# building the database-wide matrix. See compute_res_empirical_score().
C_CUTOFF = 0
E_CUTOFF = 5.0

# Loaded from build_scoring_matrix.py outputs
T_DB: int = 0
DB_FREQ: dict[str, float] = {}        # aa -> database-wide surface frequency
C_AB_MATRIX: dict[tuple, float] = {}  # (aaA, aaB) -> database-wide contact count

FALLBACK_SCORES = {
    'sc_sc': -1.0,
    'sc_mc': -0.5,
    'mc_mc': -0.5
}

CHAIN_A = "A"
CHAIN_B = "B"

ZSCORE_THRESHOLD = 2.3  # default threshold for interaction_confirmed (99% significance per 3did paper)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db_in",            required=True,
                   help="Path to input domainsplit.sqlite3")
    p.add_argument("--db_out",           required=True,
                   help="Path to output domainsplit.sqlite3")
    p.add_argument("--c_ab_matrix",           required=True,
                   help="c_ab_matrix.csv from build_scoring_matrix.py")
    p.add_argument("--db_freq",           required=True,
                   help="db_freq.csv from build_scoring_matrix.py")
    p.add_argument("--t_db",             required=True,
                   help="t_db.txt from build_scoring_matrix.py")
    # p.add_argument("--zscore_threshold", type=float, default=2.3,
    #                help="Z-score threshold for interaction_confirmed (default: 2.3, "
    #                     "i.e. 99%% significance per the 3did paper)")
    p.add_argument("--versions",         required=True,
                   help="Path to write versions.yml")
    p.add_argument("--process_name",     required=True,
                   help="Nextflow process name for versions.yml")
    p.add_argument("--source", required=True, choices=["AF3", "RF"],
                help="Select source of domain structures to use for scoring (AF3 or RF)")
    return p.parse_args()


# Load precomputed database-wide scoring artefacts from build_scoring_matrix.py outputs
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



# Check whether a single atom pair qualifies as a contact under the 3did definition
def get_atom_interaction(atomA, atomB):
    """
    Return (distance, interaction_type) for a single atom pair if it
    satisfies the 3did contact definition, else (inf, None).
    interaction_type in {'sc_sc', 'sc_mc', 'mc_mc'}.
    """
    elemA = (atomA.element or "").strip().upper()
    elemB = (atomB.element or "").strip().upper()
    pair  = tuple(sorted([elemA, elemB]))
    dist  = atomA - atomB

    qualifies = (
        (pair == ('C', 'C') and dist <= CC_VDW) or
        (pair == ('N', 'O') and dist <= NO_SB)   # covers both H-bond (<=3.5) and salt bridge (<=5.5)
    )
    if not qualifies:
        return math.inf, None

    mc_a = atomA.get_id() in BACKBONE_ATOMS
    mc_b = atomB.get_id() in BACKBONE_ATOMS
    if mc_a and mc_b:
        itype = 'mc_mc'
    elif mc_a or mc_b:
        itype = 'sc_mc'
    else:
        itype = 'sc_sc'

    return dist, itype


def check_for_interaction(resA, resB):
    """
    Check whether two residues interact under the 3did definition.
    Returns the interaction type of the closest qualifying atom pair,
    or None if no atom pair qualifies.
    """
    best_dist  = math.inf
    best_itype = None
    for atomA in resA.get_atoms():
        for atomB in resB.get_atoms():
            dist, itype = get_atom_interaction(atomA, atomB)
            if dist < best_dist:
                best_dist  = dist
                best_itype = itype
    return best_itype


def check_rsa(domain, chain_id, threshold=0.1):
    # For each residue ensure that its RSA is above the threshold (0.1) to be considered surface-exposed.
    rsa_residue_values = calculate_rsa_residue_level(domain, chain_id=chain_id)
    checked = []
    for res_id, rsa in rsa_residue_values.items():
        if rsa is not None and rsa > threshold:
            checked.append(res_id)
    return checked


def find_interacting_residues(res_a, res_b, structure):
    """
    Build the residue-pair interaction matrix for a domain pair.

    Returns:
        matrix[i][j] = interaction_type string ('sc_sc'/'sc_mc'/'mc_mc') if residues i and j interact, else []
        n_interacting = count of non-None entries
    """
    n_a, n_b = len(res_a), len(res_b)
    matrix = [[""] * n_b for _ in range(n_a)]
    # Get matrix shape
    print(len(matrix), len(matrix[0]) if matrix else 0)
    
    n_interacting = 0

    res_a_checked = check_rsa(structure, chain_id=CHAIN_A)
    res_b_checked = check_rsa(structure, chain_id=CHAIN_B)

    for i, rA in enumerate(res_a):
        if rA.get_id() not in res_a_checked:
            continue
        for j, rB in enumerate(res_b):
            if rB.get_id() not in res_b_checked:
                continue
            itype = check_for_interaction(rA, rB)
            if itype is not None:
                try:
                    matrix[i][j] = itype
                except IndexError:
                    print(f"IndexError: i={i}, j={j}, len(matrix)={len(matrix)}, len(matrix[0])={len(matrix[0]) if matrix else 0}")
                    raise
                n_interacting += 1

    return matrix, n_interacting


def get_expected_contacts(aaA, aaB):
    """
    E_ab = T_DB * freq_a * freq_b   (database-wide molar fraction random-state model).
    Computed at lookup time from the loaded T_DB / DB_FREQ rather than
    precomputed and stored, since it is a simple product of values that
    are already available after build_scoring_matrix.py has run.
    """
    return T_DB * DB_FREQ.get(aaA, 0.0) * DB_FREQ.get(aaB, 0.0)


def compute_res_empirical_score(aaA, aaB, interaction_type):
    """
    S_ab = log10(C_ab / E_ab)   if C_ab > C_CUTOFF and E_ab >= E_CUTOFF
         = fallback score        otherwise

    Application of both cutoffs
    symmetrically lookup since the database doesn't distinguish ordering.
    """
    C_ab = C_AB_MATRIX.get((aaA, aaB), C_AB_MATRIX.get((aaB, aaA), 0.0))
    E_ab = get_expected_contacts(aaA, aaB)

    if C_ab > C_CUTOFF and E_ab >= E_CUTOFF:
        return math.log10(C_ab / E_ab)
    return FALLBACK_SCORES.get(interaction_type, -1.0)


def compute_empirical_potential(residues_a, residues_b, matrix):
    """
    Sum S_ab over all interacting residue pairs (matrix[i][j] is not None).
    residues_a/b may be Bio.PDB residue objects (real complex) or plain
    3-letter AA strings (random trials) -- both are handled transparently.
    """
    total = 0.0
    for i, rA in enumerate(residues_a):
        for j, rB in enumerate(residues_b):
            itype = matrix[i][j]
            if itype is not None:
                aaA = rA if isinstance(rA, str) else rA.get_resname()
                aaB = rB if isinstance(rB, str) else rB.get_resname()
                total += compute_res_empirical_score(aaA, aaB, itype)
    return total


# Z-score setup
def sample_sequence(length, freq_dict):
    """
    Sample a random AA sequence of the given length using DATABASE-WIDE
    surface frequencies (DB_FREQ) as the sampling distribution -- not the
    frequencies of the two specific input proteins. This matches the
    "reproduce 3did scores" interpretation: the background model is the
    database, not the individual complex.
    """
    amino_acids = list(freq_dict.keys())
    weights     = list(freq_dict.values())
    total       = sum(weights)
    if total == 0:
        weights = [1.0] * len(amino_acids)
    return random.choices(amino_acids, weights=weights, k=length)


def run_random_trials(n_a, n_b, matrix, n_trials):
    """
    Run n_trials random scorings using the SAME contact geometry (matrix)
    but randomised AA identities drawn from DB_FREQ. The geometry (which
    positions interact, and at what sc/mc classification) is fixed --
    only amino acid identity is resampled each trial.
    """
    scores = []
    for _ in range(n_trials):
        rand_a = sample_sequence(n_a, DB_FREQ)
        rand_b = sample_sequence(n_b, DB_FREQ)
        scores.append(compute_empirical_potential(rand_a, rand_b, matrix))
    return scores


def compute_z_score(real_score, random_scores):
    mean_r = float(np.mean(random_scores))
    std_r  = float(np.std(random_scores))
    if std_r == 0:
        return 0.0
    return (real_score - mean_r) / std_r



def connect_db(db_path):
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA synchronous = OFF")
    conn.execute("PRAGMA journal_mode = MEMORY")
    return conn


def aggregate_scores(scores):
    """
    Aggregate score, different possible methods, return ddi_id -> confirmed
    1. Majority vote on confirmed
    2. Mean z-score -> confirmed if mean z-score >= threshold
    """
    aggregated = {}
    for ddi_id, score_list in scores.items():
        confirmed_votes = sum(confirmed for _, confirmed in score_list)
        mean_z_score = np.mean([z for z, _ in score_list if z is not None])
        # Majority vote
        majority_confirmed = int(confirmed_votes > len(score_list) / 2)
        # Mean z-score
        mean_confirmed = int(mean_z_score >= ZSCORE_THRESHOLD)
        aggregated[ddi_id] = {
            'majority_confirmed': majority_confirmed,
            'mean_confirmed': mean_confirmed
        }

    return aggregated


def update_scores_in_db(conn, aggregated_scores, source):
    """
    Update DDI table with the column interaction_confirmed based on the aggregated scores.
    NOTE: for now add for both majority and mean confirmed, but in practice we might choose one method.
    """

    new_column_names = [f'interaction_confirmed_majority_{source}', f'interaction_confirmed_mean_{source}']

    conn.execute(f"""
        ALTER TABLE domain_domain_interaction ADD COLUMN {new_column_names[0]} INTEGER DEFAULT 0;
    """)
    conn.execute(f"""
        ALTER TABLE domain_domain_interaction ADD COLUMN {new_column_names[1]} INTEGER DEFAULT 0;
    """)

    conn.executemany(f"""
        UPDATE domain_domain_interaction
        SET {new_column_names[0]} = ?,
            {new_column_names[1]} = ?
        WHERE id = ?
    """, [(scores['majority_confirmed'], scores['mean_confirmed'], ddi_id) for ddi_id, scores in aggregated_scores.items()])
    conn.commit()



# Write versions.yml for Nextflow reporting
def write_versions(path, process_name):
    import Bio
    with open(path, "w") as f:
        f.write(f'"{process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    biopython: {Bio.__version__}\n")
        f.write(f"    numpy: {np.__version__}\n")



def main():
    args = parse_args()

    print(f"[score_ddi] source={args.source}, zscore_threshold={ZSCORE_THRESHOLD}, "
          f"trials={TRIALS}", flush=True)
    

    # Load precomputed database-wide scoring artefacts
    load_c_ab_matrix(args.c_ab_matrix)
    load_db_freq(args.db_freq)
    load_t_db(args.t_db)
    print(f"  Loaded T_DB={T_DB}, {len(C_AB_MATRIX)} C_ab entries, "
          f"{len(DB_FREQ)} frequency entries", flush=True)

    # Copy DB so we never mutate the staged input
    shutil.copy(args.db_in, args.db_out)
    conn_out = connect_db(args.db_out)
    
    scores = {}

    source = args.source

    source_structures = get_domain_structures(args.db_in, source=args.source)

    print(f"[score_ddi] source={source}: {len(source_structures)} predicted "
          f"complexes to score", flush=True)


    # If no predictions exist for a source, there is nothing to do
    if len(source_structures) == 0:
        print(f"[score_ddi] WARNING: No domain structures found for source={source}, skipping scoring", flush=True)
        conn_out.close()
        write_versions(args.versions, args.process_name)
        print("[score_ddi] done", flush=True)
        return
    

        
    for (ds_id, ddi_id, pdb_gz, _source) in source_structures:

        path_to_tmp_pdb = bytes_to_tempfile(pdb_gz)
        structure = None
        try:
            structure = PDBParser(QUIET=True).get_structure(f"ddi_{ddi_id}", path_to_tmp_pdb)
        except Exception as e:
            print(f"WARNING: Failed to parse PDB for DDI {ddi_id}: {e}", file=sys.stderr)

        res_a = extract_domain_residues(structure, CHAIN_A)
        res_b = extract_domain_residues(structure, CHAIN_B)

        # Check rsa giving it a domain_structure object


        if not res_a or not res_b:
            print(f"  WARNING: no residues extracted for ddi {ddi_id}", file=sys.stderr)
            continue

        # Step 1: find interacting residue pairs and require >= 5 (3did rule)
        matrix, n_interacting = find_interacting_residues(res_a, res_b, structure)
        if n_interacting < MIN_INTERACTING_PAIRS:
            print(f"  Row {ds_id}: only {n_interacting} interacting pairs "
                    f"(< {MIN_INTERACTING_PAIRS}) -- not interacting, "
                    f"z_score=0, confirmed=0", flush=True)
            # update_score(args.db_out, ds_id, 0.0)
            continue

        # Step 2: score the real complex
        real_score = compute_empirical_potential(res_a, res_b, matrix)

        # Step 3: score random trials (same geometry, DB-wide AA sampling)
        random_scores = run_random_trials(len(res_a), len(res_b), matrix, TRIALS)

        # Step 4: z-score and significance
        z_score = compute_z_score(real_score, random_scores)
        confirmed = int(z_score >= ZSCORE_THRESHOLD)

        print(f"  Row {ds_id}: n_interacting={n_interacting}, "
                f"score={real_score:.3f}, z={z_score:.3f}, "
                f"confirmed={confirmed}", flush=True)
        # Update z-score column in domain_structure table, just to keep a record of the actual z-score for each complex
        update_score(conn_out, ds_id, z_score)
        if ddi_id not in scores:
            scores[ddi_id] = []
        scores[ddi_id].append((z_score, confirmed))
    updated_scores = aggregate_scores(scores)
    update_scores_in_db(conn_out, updated_scores, source)
    conn_out.close()


    

    write_versions(args.versions, args.process_name)
    print("[score_ddi] done", flush=True)


if __name__ == "__main__":
    main()