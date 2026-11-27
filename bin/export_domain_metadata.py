#!/usr/bin/env python3
"""Per-method domain_metadata.csv (domainbenchmark's enrichment input):
ProtParam sequence features per domain instance, ported from master_ct's
create_protein_domain_metadata.py without the InterPro counters.

Instances are those of ddi_split_membership, per method (the split databases
also hold family-context proteins, which are not benchmark instances). The
domain sequence is protein.sequence[start-1:end] (1-based, inclusive); a
window outside the sequence is fatal, since ENRICH_STRUCTURAL_AF and the
embeddings read the same coordinates.

As in the helper: a sequence with X is skipped (no row; domainbenchmark
averages per Pfam over the remaining instances), Z/B/U/O are substituted by
Q/N/C/K before computing, and average/summed flexibility is NaN for a domain
of 9 residues or fewer (ProtParam's 9-residue window yields no value). Unlike
the helper the column is spelled isoelectric_point and the sequence string is not written.

Output: <outdir>/<method>/domain_metadata.csv and --report (one TSV of counts
over the distinct instances).
"""

import argparse
import csv
import os
import sqlite3
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils_struct import extract_data_from_instance_id  # noqa: E402

SUBSTITUTIONS = {"Z": "Q", "B": "N", "U": "C", "O": "K"}
FEATURES = ["sequence_length", "molecular_weight", "isoelectric_point", "gravy_score",
            "instability_index", "aromaticity", "net_charge_ph5", "net_charge_ph7", "net_charge_ph9",
            "frac_helix", "frac_sheet", "frac_turn", "molar_ext_with_cystine", "molar_ext_with_disulfide",
            "average_flexibility", "summed_flexibility", "rel_negative_charge", "rel_positive_charge",
            "rel_no_charge", "rel_nonpolar", "rel_polar"]


def features(raw):
    """Feature dict for one domain sequence (no X)."""
    from Bio.SeqUtils import molecular_weight
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

    seq = raw.translate(str.maketrans(SUBSTITUTIONS))
    n = len(seq)
    a = ProteinAnalysis(seq)
    flex = a.flexibility() if n >= 9 else []
    helix, sheet, turn = a.secondary_structure_fraction()
    ext_reduced, ext_disulfide = a.molar_extinction_coefficient()
    frac = lambda letters: sum(seq.count(c) for c in letters) / n  # noqa: E731
    return {
        "sequence_length": n,
        "molecular_weight": molecular_weight(seq, seq_type="protein"),
        "isoelectric_point": a.isoelectric_point(),
        "gravy_score": a.gravy(),
        "instability_index": a.instability_index(),
        "aromaticity": a.aromaticity(),
        "net_charge_ph5": a.charge_at_pH(5.0),
        "net_charge_ph7": a.charge_at_pH(7.0),
        "net_charge_ph9": a.charge_at_pH(9.0),
        "frac_helix": helix,
        "frac_sheet": sheet,
        "frac_turn": turn,
        "molar_ext_with_cystine": ext_reduced,
        "molar_ext_with_disulfide": ext_disulfide,
        "average_flexibility": float(np.mean(flex)) if flex else float("nan"),
        "summed_flexibility": float(np.sum(flex)) if flex else float("nan"),
        "rel_negative_charge": frac("DE"),
        "rel_positive_charge": frac("KRH"),
        "rel_no_charge": 1.0 - frac("KRHDE"),
        "rel_nonpolar": frac("AFGILMPVW"),
        "rel_polar": frac("DEHKNQRST"),
    }


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--db", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--versions")
    p.add_argument("--process_name", default="EXPORT_DOMAIN_METADATA")
    args = p.parse_args(argv)

    conn = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    by_method = {}
    for method, ia, ib in conn.execute("SELECT method, instance_id_a, instance_id_b FROM ddi_split_membership"):
        by_method.setdefault(method, set()).update((ia, ib))
    instances = sorted(set().union(*by_method.values())) if by_method else []
    parsed = {i: extract_data_from_instance_id(i) for i in instances}
    needed = {u for _p, u, _s, _e in parsed.values()}
    sequences = {u: s for u, s in conn.execute("SELECT uniprot_id, sequence FROM protein") if u in needed}
    conn.close()

    counts = dict.fromkeys(["instances", "written", "skipped_x"] + [f"substituted_{c}" for c in SUBSTITUTIONS], 0)
    rows, bad = {}, []
    for inst in instances:
        pfam, acc, start, end = parsed[inst]
        seq = sequences.get(acc) or ""
        if not (1 <= start <= end <= len(seq)):
            bad.append(f"{inst} (sequence length {len(seq)})")
            continue
        dom = seq[start - 1:end]
        counts["instances"] += 1
        if "X" in dom:
            counts["skipped_x"] += 1
            continue
        for c in SUBSTITUTIONS:
            counts[f"substituted_{c}"] += int(c in dom)
        rows[inst] = (pfam, acc, features(dom))
        counts["written"] += 1
    if bad:
        raise SystemExit(f"{len(bad)} membership instance(s) lie outside their protein.sequence, e.g. {bad[:5]}")

    for method, members in sorted(by_method.items()):
        os.makedirs(os.path.join(args.outdir, method), exist_ok=True)
        with open(os.path.join(args.outdir, method, "domain_metadata.csv"), "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["instance_id", "pfam_id", "uniprot_id"] + FEATURES)
            for inst in sorted(members):
                if inst in rows:
                    pfam, acc, f = rows[inst]
                    w.writerow([inst, pfam, acc] + ["" if np.isnan(f[c]) else f[c] for c in FEATURES])

    with open(args.report, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        w.writerows(counts.items())
    print(f"[domain_metadata] {counts}", flush=True)

    if args.versions:
        import Bio
        with open(args.versions, "w") as f:
            f.write(f'"{args.process_name}":\n')
            f.write(f"    python: {sys.version.split()[0]}\n")
            f.write(f"    biopython: {Bio.__version__}\n")


if __name__ == "__main__":
    main()
