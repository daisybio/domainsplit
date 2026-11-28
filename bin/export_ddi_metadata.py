#!/usr/bin/env python3
"""Per-method ddi_metadata.csv (domainbenchmark's enrichment input) from the
slices ENRICH_STRUCTURAL_AF stored in structures.h5.

One row per instance pair, in membership orientation, with the 11 interface
scores of master_ct's create_ddi_metadata.py under an ``af_`` prefix. Scores
come only from the stored slice (CB coordinates, CA pLDDT) and its PAE block,
so a scoring change never re-runs slicing.

Per method, because domainbenchmark reads one metadata directory and filters
it by instance pair, while STRUCTURE_COVERAGE_BALANCE masks per (method,
split): ``<outdir>/<method>/ddi_metadata.csv`` holds a pair only if every
membership row of that pair in that method has has_structure = 1. A pair
masked in one split of a method (e.g. a positive shared by the 1:1 and 1:10
tests) is left out of the whole method. A method without any such pair gets a
header-only file, which domainbenchmark reads as "no DDI metadata".

Not a numerical reproduction of the helper:

* The slice went through a CIF -> PDB rewrite (coordinates and pLDDT rounded
  to the PDB format), so a distance at the contact cutoff can flip.
* Residues are matched to PAE rows by the residues_a / residues_b attrs that
  ENRICH_STRUCTURAL_AF wrote from the AF3 token map, not by parsed-residue
  ordinal as the helper did.
* Windows that failed ENRICH_STRUCTURAL_AF's sequence check were never sliced;
  the helper scored them at whatever residues carried those numbers.

Undefined values are NaN: a pair without any CB-CB contact has no pDockQ /
pDockQ2 / ipSAE (numberOfInteractions is 0, observed); a pair without PAE has
only min_cb_dist and pDockQ.
"""

import argparse
import csv
import gzip
import os
import sqlite3
import sys
from collections import Counter, defaultdict

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils_struct import extract_data_from_instance_id, structure_h5_keys  # noqa: E402

BASE_COLS = ["pfam_id_a", "pfam_id_b", "uniprot_id_a", "uniprot_id_b",
             "instance_id_a", "instance_id_b", "source"]
SCORE_COLS = ["min_cb_dist", "pDockQ", "pDockQ2_AB", "pDockQ2_BA", "LIS_AB", "LIS_BA",
              "ipSAE_AB", "ipTM_d0chn_AB", "ipSAE_BA", "ipTM_d0chn_BA", "numberOfInteractions"]
NOCONTACT_COLS = ["pDockQ", "pDockQ2_AB", "pDockQ2_BA", "ipSAE_AB", "ipSAE_BA"]
PAE_COLS = [c for c in SCORE_COLS if c not in ("min_cb_dist", "pDockQ")]
D0CHN = 10.0  # pDockQ2's fixed d0, as in the helper


# --------------------------------------------------------------------------
# Slice parsing and residue alignment
# --------------------------------------------------------------------------

def parse_slice(text):
    """{chain: [(resseq, ca_xyz|None, cb_xyz|None, ca_bfactor)]} in file order,
    one entry per residue, every residue kept so the list can be aligned to
    residues_a / residues_b before CA filtering. A residue is keyed like
    Biopython's residue id (hetero name, number, insertion code), so a HETATM
    sharing a number with an amino acid stays a separate entry, as it is in
    the PAE rows. Only the first alternate location seen per residue is read,
    so CA and CB never come from different conformers."""
    chains, index, altloc = defaultdict(list), {}, {}
    for line in text.splitlines():
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        het = line[17:20] if line.startswith("HETATM") else ""
        chain, key = line[21], (line[21], het, line[22:27])
        if key not in index:
            index[key] = len(chains[chain])
            chains[chain].append([int(line[22:26]), None, None, float("nan")])
        alt = line[16]
        if alt != " " and altloc.setdefault(key, alt) != alt:
            continue
        res = chains[chain][index[key]]
        atom = line[12:16].strip()
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        if atom == "CA" and res[1] is None:
            res[1], res[3] = xyz, float(line[60:66])
        elif atom == "CB" and res[2] is None:
            res[2] = xyz
    return chains


def align(chains, residues_a, residues_b, label):
    """(cb (N,3), plddt (N,), idx_a, idx_b, pae_rows (N,)) over the CA-bearing
    residues. Raises unless chain A / B list exactly residues_a / residues_b,
    in order: the PAE block's rows are [residues_a..., residues_b...]."""
    cb, plddt, idx, pae_rows = [], [], {"A": [], "B": []}, []
    offset = 0
    for chain, expected in (("A", residues_a), ("B", residues_b)):
        got = chains.get(chain, [])
        if [r[0] for r in got] != [int(x) for x in expected]:
            raise ValueError(f"{label}: chain {chain} residues do not match the PAE residue list "
                             f"({len(got)} vs {len(expected)})")
        for k, (_num, ca, cbx, bfac) in enumerate(got):
            if ca is None:  # the helper's rule: no CA, no residue
                continue
            idx[chain].append(len(cb))
            cb.append(cbx if cbx is not None else ca)
            plddt.append(bfac)
            pae_rows.append(offset + k)
        offset += len(got)
    return (np.array(cb, dtype=float).reshape(-1, 3), np.array(plddt, dtype=float),
            np.array(idx["A"], dtype=int), np.array(idx["B"], dtype=int), np.array(pae_rows, dtype=int))


def residue_numbers(chains, chain):
    return [r[0] for r in chains.get(chain, [])]


# --------------------------------------------------------------------------
# Scores (master_ct create_ddi_metadata.calculate_scores, NaN where undefined)
# --------------------------------------------------------------------------

def _d0(n):
    n = max(n, 19)
    return 1.24 * (n - 15) ** (1 / 3) - 1.8


def _pdockq2_direction(plddt, res, partner_pae):
    uniq, inv = np.unique(res, return_inverse=True)
    inv = inv.reshape(-1)
    tm = 1.0 / (1.0 + (partner_pae / D0CHN) ** 2)
    w = np.bincount(inv, weights=tm) / np.bincount(inv)
    x = np.mean(plddt[uniq] * w)
    return float(1.31034 / (1 + np.exp(-0.075 * (x - 84.733))) + 0.005)


def calculate_scores(cb, plddt, pae, idx_a, idx_b, pae_cutoff, contact_cutoff):
    """Scores for domain A (idx_a) vs domain B (idx_b). `pae` is indexed like
    cb / plddt, or None. AB = rows of A, columns of B (PAE is asymmetric)."""
    nan = float("nan")
    scores = dict.fromkeys(SCORE_COLS, nan)
    if len(idx_a) == 0 or len(idx_b) == 0:
        return scores
    diff = cb[idx_a, None, :] - cb[None, idx_b, :]
    d2 = np.einsum("ijk,ijk->ij", diff, diff)
    hi, hj = np.nonzero(d2 < contact_cutoff * contact_cutoff)
    ci, cj = idx_a[hi], idx_b[hj]
    scores["min_cb_dist"] = float(np.sqrt(d2.min()))
    d0 = _d0(float(min(len(idx_a), len(idx_b))))

    if len(ci):
        interface = np.unique(np.concatenate([ci, cj]))
        x = plddt[interface].mean() * np.log10(len(ci))
        scores["pDockQ"] = float(0.724 / (1 + np.exp(-0.052 * (x - 152.611))) + 0.018)
    if pae is None:
        return scores

    for direction, sub in (("AB", pae[np.ix_(idx_a, idx_b)]), ("BA", pae[np.ix_(idx_b, idx_a)])):
        mask = sub < pae_cutoff
        scores[f"LIS_{direction}"] = float((1.0 - sub[mask] / pae_cutoff).mean()) if mask.any() else 0.0
        scores[f"ipTM_d0chn_{direction}"] = float((1.0 / (1.0 + (sub / d0) ** 2)).max(axis=1).mean())
    if len(ci) == 0:
        scores["numberOfInteractions"] = 0
        return scores
    pae_ab, pae_ba = pae[ci, cj], pae[cj, ci]
    scores["pDockQ2_AB"] = _pdockq2_direction(plddt, ci, pae_ab)
    scores["pDockQ2_BA"] = _pdockq2_direction(plddt, cj, pae_ba)
    for direction, vals in (("AB", pae_ab), ("BA", pae_ba)):
        passing = vals[vals < pae_cutoff]
        scores[f"ipSAE_{direction}"] = float(np.mean(1.0 / (1.0 + (passing / d0) ** 2))) if len(passing) else 0.0
    scores["numberOfInteractions"] = int(np.sum((pae_ab < pae_cutoff) | (pae_ba < pae_cutoff)))
    return scores


def score_pair(h5, ia, ib, pae_cutoff, contact_cutoff):
    """(source, scores, has_pae) for one stored slice."""
    group, dset = structure_h5_keys(ia, ib)
    node = h5[group][dset]
    source = node.attrs.get("af_source", "")
    source = source.decode() if isinstance(source, bytes) else str(source)
    chains = parse_slice(gzip.decompress(node[()].tobytes()).decode("utf-8"))
    pae_node = h5.get(f"pae/{group}/{dset}")
    label = f"{group}/{dset}"
    if pae_node is None:
        res_a, res_b = residue_numbers(chains, "A"), residue_numbers(chains, "B")
    else:
        res_a, res_b = pae_node.attrs["residues_a"], pae_node.attrs["residues_b"]
    cb, plddt, idx_a, idx_b, rows = align(chains, res_a, res_b, label)
    pae = None
    if pae_node is not None:
        full = np.asarray(pae_node[()], dtype=float)
        n = len(res_a) + len(res_b)
        if full.shape != (n, n):
            raise ValueError(f"{label}: PAE block {full.shape} does not fit {n} residues")
        pae = full[np.ix_(rows, rows)]
    return source, calculate_scores(cb, plddt, pae, idx_a, idx_b, pae_cutoff, contact_cutoff), pae is not None


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def membership(db):
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    rows = conn.execute(
        "SELECT m.method, m.split, m.instance_id_a, m.instance_id_b, d.negative, "
        "m.has_structure, m.has_structure_raw "
        "FROM ddi_split_membership m JOIN domain_domain_interaction d ON d.id = m.ddi_id").fetchall()
    conn.close()
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--db", required=True, help="master after STRUCTURE_COVERAGE_BALANCE")
    p.add_argument("--structures_h5", default="", help="ENRICH_STRUCTURAL_AF's structures.h5; empty = no structure")
    p.add_argument("--pae_cutoff", type=float, default=12.0)
    p.add_argument("--contact_cutoff", type=float, default=8.0)
    p.add_argument("--outdir", required=True)
    p.add_argument("--coverage", required=True)
    p.add_argument("--versions")
    p.add_argument("--process_name", default="EXPORT_DDI_METADATA")
    args = p.parse_args(argv)

    rows = membership(args.db)
    # Per (method, pair): every row structured -> exportable for that method.
    all_structured = {}
    for method, _split, ia, ib, _neg, has, _raw in rows:
        key = (method, ia, ib)
        all_structured[key] = all_structured.get(key, True) and bool(has)
    exported = {k for k, ok in all_structured.items() if ok}
    pairs = sorted({(ia, ib) for _m, ia, ib in exported})
    if pairs and not args.structures_h5:
        raise SystemExit(f"{len(pairs)} instance pair(s) have has_structure = 1 but no --structures_h5 was given")

    scored = {}  # (ia, ib) -> (source, scores, has_pae)
    if args.structures_h5:
        import h5py
        with h5py.File(args.structures_h5, "r") as h5:
            for n, (ia, ib) in enumerate(pairs, 1):
                group, dset = structure_h5_keys(ia, ib)
                if group not in h5 or dset not in h5[group]:
                    raise SystemExit(f"{ia}/{ib}: has_structure = 1 but no slice in {args.structures_h5}")
                scored[(ia, ib)] = score_pair(h5, ia, ib, args.pae_cutoff, args.contact_cutoff)
                if n % 5000 == 0:
                    print(f"[ddi_metadata] scored {n}/{len(pairs)} instance pairs", flush=True)

    methods = sorted({r[0] for r in rows})
    for method in methods:
        out_dir = os.path.join(args.outdir, method)
        os.makedirs(out_dir, exist_ok=True)
        with open(os.path.join(out_dir, "ddi_metadata.csv"), "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(BASE_COLS + [f"af_{c}" for c in SCORE_COLS])
            for ia, ib in pairs:
                if (method, ia, ib) not in exported:
                    continue
                source, scores, _ = scored[(ia, ib)]
                pa, ua, _, _ = extract_data_from_instance_id(ia)
                pb, ub, _, _ = extract_data_from_instance_id(ib)
                # NaN as an empty field: pandas reads it back as NaN.
                w.writerow([pa, pb, ua, ub, ia, ib, source]
                           + ["" if np.isnan(scores[c]) else scores[c] for c in SCORE_COLS])

    # Coverage per membership row: method x split x label x source.
    cov = Counter()
    for method, split, ia, ib, neg, has, raw in rows:
        label = "negative" if neg else "positive"
        hit = scored.get((ia, ib)) if (method, ia, ib) in exported else None
        source = hit[0] if hit else "none"
        k = (method, split, label, source)
        cov[k + ("rows",)] += 1
        cov[k + ("with_structure_raw",)] += int(bool(raw))
        cov[k + ("with_structure",)] += int(bool(has))
        if hit:
            cov[k + ("exported",)] += 1
            cov[k + ("with_pae",)] += int(hit[2])
            cov[k + ("with_contact",)] += int(not np.isnan(hit[1]["pDockQ"]))
    fields = ["rows", "with_structure_raw", "with_structure", "exported", "with_pae", "with_contact"]
    with open(args.coverage, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["method", "split", "label", "source"] + fields)
        for k in sorted({k[:4] for k in cov}):
            w.writerow(list(k) + [cov[k + (f,)] for f in fields])
    print(f"[ddi_metadata] {len(scored)} instance pairs scored, {len(methods)} method file(s) written", flush=True)

    if args.versions:
        with open(args.versions, "w") as f:
            f.write(f'"{args.process_name}":\n')
            f.write(f"    python: {sys.version.split()[0]}\n")
            f.write(f"    numpy: {np.__version__}\n")


if __name__ == "__main__":
    main()
