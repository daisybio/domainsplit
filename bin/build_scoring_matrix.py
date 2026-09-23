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

SHARDING: when --n_shards > 1, this script processes only its slice of
the method's train structures and writes RAW, UN-normalized counts:
    <method>.<shard_id>.c_ab_counts.csv     — 20x20 contact counts (summable)
    <method>.<shard_id>.surface_counts.csv  — surface AA counts (summable)
This is a reduction (sums), so any split across shards is safe -- unlike
score_ddi.py, no grouping is needed here. Run merge_scoring_matrix.py
once all shards for a method are done to sum these and produce the final,
normalized:
    <method>.c_ab_matrix.csv — 20x20 contact counts (aaA, aaB, count)
    <method>.db_freq.csv     — train-set surface AA frequencies (aa, frequency)
    <method>.t_db.txt        — total contact count T_DB
which score_ddi.py consumes unchanged. With --n_shards=1 (default), a
single "merge" of one shard still runs, so behaviour/output is identical
to before, just via the two-step shard+merge path.

PERFORMANCE NOTE (vs. original):
  - Contact detection (residues_contact) is vectorized with NumPy/cdist
    instead of a Python double-loop over every atom pair per residue pair.
    Only existence of a qualifying atom pair is needed here (not the
    closest one, or its mc/sc type, as in score_ddi.py), so this is a
    lighter version of the same fix.
  - PDB parsing reads structures from an in-memory buffer instead of a
    temp file on disk (see utils_struct.bytes_to_stringio).
"""
import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import Bio
import h5py
import numpy as np
import scipy
from scipy.spatial.distance import cdist
from Bio.PDB.PDBParser import PDBParser

from utils_struct import (
    AA_3, bytes_to_stringio, calculate_rsa_residue_level,
    extract_domain_residues, get_positive_train_structures,
    get_structure_bytes,
)

CC_VDW = 5.0   # vdW:        C-C distance <= 5.0 Å
NO_SB  = 5.5   # Salt bridge: N-O distance <= 5.5 Å

CHAIN_A, CHAIN_B = "A", "B"
# Relative solvent accessibility of the unbound protein must be >= 10% to
# count a residue as surface (excludes buried side-chains, per Sprinzak &
# Margalit 2001). Applied to both the contact matrix and DB_FREQ, since
# both are meant to reflect surface, not buried, residues.
RSA_THRESHOLD = 0.1


def check_rsa(domain, chain_id, threshold=RSA_THRESHOLD):
    """Residue ids with relative solvent accessibility >= threshold."""
    rsa = calculate_rsa_residue_level(domain, chain_id=chain_id)
    return {rid for rid, v in rsa.items() if v is not None and v >= threshold}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db_in", required=True)
    p.add_argument("--structures_h5", required=True)
    p.add_argument("--method", required=True)
    p.add_argument("--c_ab_counts_out", required=True,
                    help="raw (un-normalized) 20x20 contact-count output for this shard")
    p.add_argument("--surface_counts_out", required=True,
                    help="raw (un-normalized) surface AA count output for this shard")
    p.add_argument("--versions", required=True)
    p.add_argument("--process_name", required=True)
    p.add_argument("--shard_id", type=int, default=0,
                    help="0-based index of this shard (default 0 = no sharding)")
    p.add_argument("--n_shards", type=int, default=1,
                    help="total number of shards this method's train structures are split into")
    args = p.parse_args()
    if not (0 <= args.shard_id < args.n_shards):
        p.error(f"--shard_id ({args.shard_id}) must be in [0, --n_shards={args.n_shards})")
    return args


def _residue_atom_arrays(residues, checked_ids):
    """
    Flatten C/N/O atoms of the given residues (restricted to those in
    checked_ids, i.e. RSA-passing / surface residues) into parallel NumPy
    arrays:
        coords  (n_atoms, 3) float
        elems   (n_atoms,)   'C'/'N'/'O'
        res_idx (n_atoms,)   int index into `residues`
    Only C/N/O atoms are kept since those are the only elements the
    contact definition (C-C vdW, N-O salt bridge) ever uses.
    """
    coords, elems, res_idx = [], [], []
    for i, r in enumerate(residues):
        if r.get_id() not in checked_ids:
            continue
        for atom in r.get_atoms():
            el = (atom.element or "").strip().upper()
            if el not in ("C", "N", "O"):
                continue
            coords.append(atom.coord)
            elems.append(el)
            res_idx.append(i)
    return (
        np.asarray(coords, dtype=float).reshape(-1, 3),
        np.asarray(elems),
        np.asarray(res_idx, dtype=int),
    )


def find_contact_pairs_vectorized(res_a, res_b, checked_a, checked_b):
    """
    Vectorized replacement for the residues_contact() double-loop.
    Returns the set of unique (i, j) residue-index pairs (into res_a,
    res_b respectively) that have at least one qualifying atom contact,
    matching the original 3did contact definition (C-C vdW or N-O salt
    bridge), restricted to residues with RSA >= RSA_THRESHOLD on each
    side. Distance/order beyond "does a qualifying contact exist" is
    irrelevant here, unlike in score_ddi.py.
    """
    coords_a, elems_a, ridx_a = _residue_atom_arrays(res_a, checked_a)
    coords_b, elems_b, ridx_b = _residue_atom_arrays(res_b, checked_b)

    if len(coords_a) == 0 or len(coords_b) == 0:
        return set()

    dist = cdist(coords_a, coords_b)

    cc_mask = (elems_a[:, None] == "C") & (elems_b[None, :] == "C") & (dist <= CC_VDW)
    no_mask = (
        ((elems_a[:, None] == "N") & (elems_b[None, :] == "O")) |
        ((elems_a[:, None] == "O") & (elems_b[None, :] == "N"))
    ) & (dist <= NO_SB)
    hit = cc_mask | no_mask
    if not hit.any():
        return set()

    ai, bi = np.where(hit)
    ra = ridx_a[ai]
    rb = ridx_b[bi]
    return set(zip(ra.tolist(), rb.tolist()))


def build_matrix(args):
    C_ab = defaultdict(lambda: defaultdict(int))
    surface_counts = defaultdict(int)
    n_instances = n_skipped = n_no_contact = n_missing_bytes = 0

    train_structures = get_positive_train_structures(args.db_in, args.method)
    print(f"[build_scoring_matrix] method={args.method}: {len(train_structures)} "
          f"train-split positive-DDI structures total", flush=True)

    if args.n_shards > 1:
        # A pure reduction over independent structures -- any split is
        # safe to sum back together later, unlike score_ddi.py's per-ddi_id
        # aggregation. Simple round-robin by structure index is fine here.
        train_structures = train_structures[args.shard_id::args.n_shards]
        print(f"[build_scoring_matrix] method={args.method}: shard {args.shard_id}/{args.n_shards} "
              f"-> {len(train_structures)} structures", flush=True)

    with h5py.File(args.structures_h5, "r") as h5file:
        for (ddi_id, instance_id_a, instance_id_b) in train_structures:
            n_instances += 1
            pdb_gz = get_structure_bytes(h5file, instance_id_a, instance_id_b)
            if pdb_gz is None:
                print(f"WARNING: no structure bytes for ddi {ddi_id} "
                      f"({instance_id_a}, {instance_id_b})", file=sys.stderr)
                n_missing_bytes += 1
                continue

            structure = None
            try:
                structure = PDBParser(QUIET=True).get_structure(f"ddi_{ddi_id}", bytes_to_stringio(pdb_gz))
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

            # RSA of the unbound domain must be >= 10% to count as surface
            # (excludes buried side-chains); applied to both the surface
            # frequency counts and the contact matrix.
            checked_a = check_rsa(structure, chain_id=CHAIN_A)
            checked_b = check_rsa(structure, chain_id=CHAIN_B)

            for r in res_a:
                if r.get_id() in checked_a:
                    surface_counts[r.get_resname()] += 1
            for r in res_b:
                if r.get_id() in checked_b:
                    surface_counts[r.get_resname()] += 1

            contact_pairs = find_contact_pairs_vectorized(res_a, res_b, checked_a, checked_b)
            for i, j in contact_pairs:
                aaA = res_a[i].get_resname()
                aaB = res_b[j].get_resname()
                C_ab[aaA][aaB] += 1
                C_ab[aaB][aaA] += 1
            if not contact_pairs:
                n_no_contact += 1
            if n_instances % 500 == 0:
                print(f"  Processed {n_instances} DDIs...", flush=True)

    print(f"Done. Instances: {n_instances}, skipped: {n_skipped}, "
          f"missing structure bytes: {n_missing_bytes}, no contacts: {n_no_contact}", flush=True)
    return C_ab, surface_counts


def write_c_ab_counts(C_ab, path):
    """Raw (un-normalized) contact counts for this shard -- summable across
    shards by merge_scoring_matrix.py. Same format as the old final
    c_ab_matrix.csv, since it was already raw counts, not frequencies."""
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["aaA", "aaB", "count"])
        for aaA in AA_3:
            for aaB in AA_3:
                writer.writerow([aaA, aaB, C_ab[aaA][aaB]])
    print(f"Written raw C_ab counts to {path}", flush=True)


def write_surface_counts(surface_counts, path):
    """Raw (un-normalized) surface AA counts for this shard -- summable
    across shards. merge_scoring_matrix.py normalizes into db_freq.csv."""
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["aa", "count"])
        for aa in AA_3:
            writer.writerow([aa, surface_counts.get(aa, 0)])
    print(f"Written raw surface counts to {path}", flush=True)


def write_versions(path, process_name):
    with open(path, "w") as f:
        f.write(f'"{process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    biopython: {Bio.__version__}\n")
        f.write(f"    numpy: {np.__version__}\n")
        f.write(f"    scipy: {scipy.__version__}\n")
        f.write(f"    h5py: {h5py.__version__}\n")


def main():
    args = parse_args()
    print(f"[build_scoring_matrix] method={args.method}: "
          f"building shard {args.shard_id}/{args.n_shards} of the train-set scoring matrix", flush=True)
    C_ab, surface_counts = build_matrix(args)
    write_c_ab_counts(C_ab, args.c_ab_counts_out)
    write_surface_counts(surface_counts, args.surface_counts_out)
    write_versions(args.versions, args.process_name)
    print(f"[build_scoring_matrix] method={args.method}: shard {args.shard_id}/{args.n_shards} done", flush=True)


if __name__ == "__main__":
    main()
