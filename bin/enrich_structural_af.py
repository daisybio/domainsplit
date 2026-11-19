#!/usr/bin/env python3
"""Slice every AF3-predicted complex needed by ddi_split_membership into its
domain-instance-pair structures, store them in structures.h5, and
set ddi_split_membership.has_structure = 1 on every row of each sliced
(ddi_id, instance_id_a, instance_id_b) triple.

Driven by ddi_split_membership: only instance pairs that one of --methods
uses in some split are sliced. Grouped by UNORDERED protein pair so a complex predicted once
is loaded/parsed exactly once, even when it realizes several DDIs and even
when the pair shows up as (A, B) for some DDIs and (B, A) for others.

When af_metadata has no entry for a pair, or the referenced file is missing
or unreadable, nothing is written for it and has_structure keeps its default
of 0. The same holds for an instance pair whose window does not match this
run's UniProt sequence (protein.sequence) by residue number: af_metadata is
user-generated, and a complex predicted from an isoform, an older sequence
version or a renumbered construct would otherwise be sliced at the wrong
residues (see check_window). has_structure is a property of the instance pair, so it is set for
every method's rows of that triple, not only the --methods ones: the
structure is in structures.h5 either way.

AF confidence is stored next to each slice (see load_confidence): scalar
attrs iptm/ptm/ranking_score/af_source[/chain_pair_iptm] on the structure
dataset, and the PAE submatrix over the two sliced windows (rows/cols =
[domain A residues, domain B residues], row orientation) at
pae/<group_key>/<dset_key>. A missing confidence file leaves NaN / no PAE
and is counted; a present but unparsable or inconsistent one is fatal.

Every --methods membership row without structure gets a reason, reported in
--coverage (per method, split, label) and --af3_to_model_dir (one TSV per
method plus all_methods.tsv: the protein pairs to run AF3 on).
"""
import argparse as ap
import gzip
import io
import itertools
import json
import multiprocessing
import os
import re
import shutil
import sqlite3
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

import Bio
import h5py
import numpy as np
import pandas as pd
import zstandard
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.SeqUtils import seq1

import utils_struct

# How many (ddi_id, instance_id_a, instance_id_b) rows to batch per
# executemany() + commit. Keeps a single Python list from growing unbounded on
# a very large run while still avoiding one round-trip/commit per row.
INSERT_BATCH_SIZE = 20000

# Columns handed to the workers, in this order.
ROW_COLS = ["ddi_id", "instance_id_a", "instance_id_b",
            "uniprot_id_a", "start_a", "end_a",
            "uniprot_id_b", "start_b", "end_b"]

# A window must have at least this fraction of its residues present in the
# complex AND matching the UniProt sequence at their number: a truncated
# construct or a model missing most of the domain is not this domain.
MIN_WINDOW_COVERAGE = 0.8


def default_cpus() -> int:
    try:
        return len(os.sched_getaffinity(0))  # respects SLURM/cgroup pinning
    except AttributeError:
        return os.cpu_count() or 1


def parse_args():
    p = ap.ArgumentParser(description=__doc__, formatter_class=ap.RawDescriptionHelpFormatter)
    p.add_argument("--db_in", required=True)
    p.add_argument("--db_out", required=True)
    p.add_argument("--structures_h5", required=True)
    p.add_argument("--af_metadata", required=True)
    p.add_argument("--methods", required=True,
                   help="comma-separated methods whose instance pairs are sliced")
    p.add_argument("--max_mismatch_rate", type=float, required=True,
                   help="max fraction of a window's residues whose type differs from "
                        "protein.sequence at their number (as afdb_max_mismatch_rate)")
    p.add_argument("--coverage", required=True,
                   help="TSV: structure coverage per (method, split, label)")
    p.add_argument("--af3_to_model_dir", required=True,
                   help="directory for <method>.tsv + all_methods.tsv of protein pairs without structure")
    p.add_argument("--versions", required=True)
    p.add_argument("--process_name", required=True)
    p.add_argument("--cpus", type=int, default=None,
                   help="Number of worker processes (default: all CPUs available to this job). "
                        "1 = run everything in the main process (handy for debugging).")
    return p.parse_args()


def build_meta_index(meta_path: str) -> dict:
    """key: frozenset({uniprot_id_a, uniprot_id_b})
    value: (path, file_uniprot_a, source) -- file_uniprot_a is whichever
    protein is chain A in the predicted file (af_metadata's own
    uniprot_id_a), so callers can tell whether their own (a, b) matches the
    file's orientation or is swapped."""
    meta = pd.read_csv(meta_path, sep=",", dtype=str)
    # A relative path is relative to the CSV's real location: Nextflow stages
    # the CSV as a symlink, and a fixture cannot know the absolute path.
    base = os.path.dirname(os.path.realpath(meta_path))
    # Row numbers as a spreadsheet shows them: header is row 1.
    blank = meta.index[meta["path"].isna() | (meta["path"].str.strip() == "")]
    if len(blank):
        sys.exit(f"af_metadata {meta_path}: blank path on row(s) "
                 f"{', '.join(str(i + 2) for i in blank)}")
    index, first_row = {}, {}
    for i, row in zip(meta.index, meta.itertuples(index=False)):
        key = frozenset({row.uniprot_id_a, row.uniprot_id_b})
        entry = (os.path.join(base, row.path), row.uniprot_id_a, row.source)
        if key in index and index[key] != entry:
            # No rule says which of two predictions of one pair is right, so neither wins.
            sys.exit(f"af_metadata {meta_path}: rows {first_row[key]} and {i + 2} give "
                     f"conflicting (path, chain-A protein, source) for protein pair "
                     f"{'/'.join(sorted(key))}")
        index.setdefault(key, entry)
        first_row.setdefault(key, i + 2)
    return index


# ---------------------------------------------------------------------------
# Worker-side code (runs in child processes)
# ---------------------------------------------------------------------------

def read_zst(path: str) -> bytes:
    """Stream-decompress a .zst file (frames need not record their content size)."""
    chunks = []
    with open(path, "rb") as fh:
        reader = zstandard.ZstdDecompressor().stream_reader(fh)
        while True:
            chunk = reader.read(1 << 20)
            if not chunk:
                break
            chunks.append(chunk)
    return b"".join(chunks)


def load_structure(path: str, source: str):
    """Return (Biopython Structure, PDB text or None), fully in memory (no scratch files).
       source 'downloaded': `path` is a .pdb.zst AFDB model, decompressed on the fly
       source 'inferred':   `path` is a .cif file"""
    if source == "downloaded":
        text = read_zst(path).decode("utf-8")
        return PDBParser(QUIET=True).get_structure("complex", io.StringIO(text)), text
    if source == "inferred":
        return MMCIFParser(QUIET=True).get_structure("complex", path), None
    raise ValueError(f"Unknown source {source} for path {path}")


# ---------------------------------------------------------------------------
# AF confidence (ipTM / pTM / ranking_score / PAE)
# ---------------------------------------------------------------------------

AF3_MODEL_SUFFIX = "_model.cif"
AFDB_MODEL_RE = re.compile(r"-model_v(\d+)\.pdb\.zst$")
# shortcut: format not verified against a real AFDB heterodimer PDB; an
# ipTM in a REMARK line is read if present, otherwise ipTM stays NaN.
HEADER_IPTM_RE = re.compile(r"\bip_?tm\b\D{0,5}?(\d*\.\d+|\d+)", re.IGNORECASE)


def confidence_paths(path: str, source: str):
    """(summary_json, pae_json) siblings of a model file; None where the
    layout has no such file or `path` is not in the expected form."""
    if source == "inferred" and path.endswith(AF3_MODEL_SUFFIX):
        stem = path[:-len(AF3_MODEL_SUFFIX)]
        return stem + "_summary_confidences.json", stem + "_confidences.json"
    if source == "downloaded" and AFDB_MODEL_RE.search(path):
        return None, AFDB_MODEL_RE.sub(r"-predicted_aligned_error_v\1.json.zst", path)
    return None, None


def read_json(path: str):
    try:
        raw = read_zst(path) if path.endswith(".zst") else open(path, "rb").read()
        return json.loads(raw)
    except Exception as exc:
        raise ValueError(f"{path}: unparsable confidence file: {exc!r}") from exc


def _float_or_nan(value, what):
    if value is None:
        return float("nan")
    try:
        return float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{what}: not a number: {value!r}")


def header_iptm(text):
    for line in (text or "").splitlines():
        if line.startswith("REMARK"):
            m = HEADER_IPTM_RE.search(line)
            if m:
                return float(m.group(1))
    return float("nan")


def load_confidence(path: str, source: str, model, pdb_text):
    """Per protein pair: {"scalars", "chain_pair_iptm", "pae", "index", "missing"}.
    index maps (chain_id, residue.id) of `model` to a PAE row (None = no unique
    token). Missing files land in "missing"; a present file that does not
    parse or does not fit `model` raises ValueError (fatal for the run)."""
    summary_path, pae_path = confidence_paths(path, source)
    out = {"scalars": {"iptm": float("nan"), "ptm": float("nan"), "ranking_score": float("nan")},
           "chain_pair_iptm": None, "pae": None, "index": None, "missing": []}
    if source == "downloaded":
        out["scalars"]["iptm"] = header_iptm(pdb_text)
    else:
        if summary_path is None or not os.path.exists(summary_path):
            out["missing"].append(summary_path or f"{path} (summary_confidences not derivable)")
        else:
            summary = read_json(summary_path)
            if not isinstance(summary, dict):
                raise ValueError(f"{summary_path}: expected a JSON object")
            for k in out["scalars"]:
                out["scalars"][k] = _float_or_nan(summary.get(k), f"{summary_path} {k}")
            cpi = summary.get("chain_pair_iptm")
            if cpi is not None:
                try:
                    cpi = np.array(cpi, dtype=float)
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"{summary_path}: chain_pair_iptm not numeric: {exc!r}")
                if cpi.shape != (2, 2) or not np.isfinite(cpi).all():
                    raise ValueError(f"{summary_path}: chain_pair_iptm is not a finite 2x2 matrix: {cpi.tolist()}")
                out["chain_pair_iptm"] = cpi

    if pae_path is None or not os.path.exists(pae_path):
        out["missing"].append(pae_path or f"{path} (PAE file not derivable)")
        return out
    data = read_json(pae_path)
    residues = [(chain.id, r.id) for chain in model for r in chain]
    try:
        if source == "inferred":
            pae = np.array(data["pae"], dtype=float)
            chain_ids, res_ids = list(data["token_chain_ids"]), [int(x) for x in data["token_res_ids"]]
            if len(chain_ids) != len(res_ids):
                raise ValueError(f"{len(chain_ids)} token_chain_ids vs {len(res_ids)} token_res_ids")
            tokens = {}
            for i, key in enumerate(zip(chain_ids, res_ids)):
                tokens.setdefault(key, []).append(i)
            # A modified residue / ligand spans several tokens: no single PAE row.
            index = {(c, rid): (tokens[(c, rid[1])][0]
                                if rid[2] == " " and len(tokens.get((c, rid[1]), ())) == 1 else None)
                     for c, rid in residues}
            n = len(chain_ids)
        else:
            entry = data[0] if isinstance(data, list) else data
            pae = np.array(entry["predicted_aligned_error"], dtype=float)
            # AFDB PAE rows follow the PDB's chain order and residue order.
            index = {key: i for i, key in enumerate(residues)}
            n = len(residues)
    except (KeyError, IndexError, TypeError, ValueError) as exc:
        raise ValueError(f"{pae_path}: unparsable confidence file: {exc!r}") from exc
    if pae.shape != (n, n):
        raise ValueError(f"{pae_path}: PAE shape {pae.shape} does not fit the model's {n} "
                         f"{'tokens' if source == 'inferred' else 'residues'} ({path})")
    # A JSON null becomes NaN under dtype=float and would pass the shape check.
    if not np.isfinite(pae).all():
        raise ValueError(f"{pae_path}: PAE has non-finite entries ({path})")
    out["pae"], out["index"] = pae, index
    return out


def slice_confidence(conf, model, chain_a, start_a, end_a, chain_b, start_b, end_b, source):
    """(attrs, pae_sub or None, residues_a, residues_b) for one slice, in the
    row's orientation: chain_a is the row's protein A whatever its file chain."""
    attrs = dict(conf["scalars"], af_source=source)
    if conf["chain_pair_iptm"] is not None:
        order = [0, 1] if chain_a == "A" else [1, 0]
        attrs["chain_pair_iptm"] = conf["chain_pair_iptm"][np.ix_(order, order)]
    if conf["pae"] is None:
        return attrs, None, None, None
    rows, resnums = [], []
    for chain, start, end in ((chain_a, start_a, end_a), (chain_b, start_b, end_b)):
        sel = [r for r in model[chain] if start <= r.id[1] <= end]  # what DdiPairSelect keeps
        idx = [conf["index"].get((chain, r.id)) for r in sel]
        if any(i is None for i in idx):
            raise ValueError(f"chain {chain} {start}-{end}: residue(s) without a unique PAE token")
        rows += idx
        resnums.append(np.array([r.id[1] for r in sel], dtype=np.int32))
    pae = conf["pae"][np.ix_(rows, rows)].astype(np.float16)
    return attrs, pae, resnums[0], resnums[1]


def check_window(chain, start: int, end: int, sequence, max_mismatch_rate: float):
    """Raise ValueError unless `chain`'s residues numbered start..end are
    `sequence`'s residues at those numbers. Compared by residue number, not
    by order, so a gap in the model does not shift the comparison. A residue
    numbered past the sequence end counts as a mismatch. An insertion-coded
    residue (10A) has no UniProt position, so it is never validated."""
    if not sequence:
        raise ValueError("no canonical sequence to validate against")
    observed = [r for r in chain if r.id[0] == " " and r.id[2] == " " and start <= r.id[1] <= end]
    mismatches = sum(1 for r in observed
                     if r.id[1] > len(sequence) or seq1(r.get_resname()) != sequence[r.id[1] - 1])
    if observed and mismatches / len(observed) > max_mismatch_rate:
        raise ValueError(f"chain {chain.id} {start}-{end}: {mismatches}/{len(observed)} residues "
                         f"differ from the UniProt sequence (> {max_mismatch_rate})")
    coverage = (len(observed) - mismatches) / (end - start + 1)
    if coverage < MIN_WINDOW_COVERAGE:
        raise ValueError(f"chain {chain.id} {start}-{end}: only {coverage:.0%} of the window "
                         f"validated (< {MIN_WINDOW_COVERAGE:.0%})")


def slice_pair(structure, chain_a: str, start_a: int, end_a: int,
               chain_b: str, start_b: int, end_b: int) -> bytes:
    """Slicing a pair of domains from a structure. The domain from `chain_a` is written
    as chain A and the one from `chain_b` as chain B; the chain ids are
    restored afterwards so the structure can be reused for the next slice."""
    model = structure[0]
    ca, cb = model[chain_a], model[chain_b]
    # Two-step rename via temporary IDs: Biopython's Entity.id setter raises if
    # the new id already exists in the parent, so a direct A<->B swap would fail.
    ca.id, cb.id = "X", "Y"
    ca.id, cb.id = "A", "B"
    try:
        # PDBIO happily writes an empty slice; a pair needs residues on both sides.
        for chain, start, end in (("A", start_a, end_a), ("B", start_b, end_b)):
            if not any(start <= res.id[1] <= end for res in model[chain]):
                raise ValueError(f"no residues of chain {chain} in {start}-{end}")
        buf = io.StringIO()
        io_obj = PDBIO()
        io_obj.set_structure(structure)
        io_obj.save(buf, utils_struct.DdiPairSelect("A", start_a, end_a, "B", start_b, end_b))
        return gzip.compress(buf.getvalue().encode("utf-8"))
    finally:
        ca.id, cb.id = "X", "Y"
        ca.id, cb.id = chain_a, chain_b


def process_pair(task):
    """One protein pair: load the complex once, slice every instance pair."""
    path, source, file_uniprot_a, rows, sequences, max_mismatch_rate = task
    try:
        structure, pdb_text = load_structure(path, source)
    except Exception as exc:  # missing / corrupt file -> no structure
        return {"unreadable_rows": [(r[0], r[1], r[2]) for r in rows],
                "unreadable_reason": "file_missing" if not os.path.exists(path) else "file_unreadable",
                "error": f"{path}: {exc!r}",
                "slices": {}, "confidence": {}, "conf_missing": [], "ok_rows": [], "failed": []}
    # Outside any try: an unparsable confidence file must stop the run.
    conf = load_confidence(path, source, structure[0], pdb_text)

    slices = {}    # (group_key, dset_key) -> gzip PDB bytes
    confidence = {}  # (group_key, dset_key) -> slice_confidence() tuple
    errors = {}    # (group_key, dset_key) -> error string
    ok_rows, failed = [], []
    for (ddi_id, ia, ib, uni_a, sa, ea, uni_b, sb, eb) in rows:
        key = utils_struct.structure_h5_keys(ia, ib)
        if key not in slices and key not in errors:
            # Row's protein A is chain A in the file unless the file has the pair swapped.
            chain_a, chain_b = ("B", "A") if file_uniprot_a != uni_a else ("A", "B")
            model = structure[0]
            try:
                check_window(model[chain_a], int(sa), int(ea), sequences.get(uni_a), max_mismatch_rate)
                check_window(model[chain_b], int(sb), int(eb), sequences.get(uni_b), max_mismatch_rate)
                slices[key] = slice_pair(structure, chain_a, int(sa), int(ea),
                                         chain_b, int(sb), int(eb))
            except Exception as exc:
                errors[key] = f"{exc!r} (from {path})"
            else:
                try:
                    confidence[key] = slice_confidence(conf, model, chain_a, int(sa), int(ea),
                                                       chain_b, int(sb), int(eb), source)
                except ValueError as exc:
                    raise ValueError(f"{path} {ia}/{ib}: {exc}") from exc
        if key in errors:
            failed.append((ddi_id, ia, ib, errors[key]))
        else:
            ok_rows.append((ddi_id, ia, ib))
    return {"unreadable_rows": [], "unreadable_reason": None, "error": None,
            "slices": slices, "confidence": confidence, "conf_missing": conf["missing"],
            "ok_rows": ok_rows, "failed": failed}


# ---------------------------------------------------------------------------
# Main-process code
# ---------------------------------------------------------------------------

def run_tasks(tasks, n_workers: int):
    """Yield process_pair results as they finish. With n_workers > 1 at most
    4*n_workers tasks are in flight, so results can't pile up in memory faster
    than the main process writes them out."""
    if n_workers <= 1:
        for t in tasks:
            yield process_pair(t)
        return

    # 'spawn' so children never inherit the open HDF5/sqlite handles.
    ctx = multiprocessing.get_context("spawn")
    it = iter(tasks)
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as ex:
        pending = {ex.submit(process_pair, t) for t in itertools.islice(it, 4 * n_workers)}
        while pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for fut in done:
                yield fut.result()
                nxt = next(it, None)
                if nxt is not None:
                    pending.add(ex.submit(process_pair, nxt))


def fetch_relevant_instances(conn, methods) -> pd.DataFrame:
    """Every (ddi_id, instance_id_a, instance_id_b) triple used by `methods`
    in any split, deduplicated, with pfam/uniprot/
    coordinates parsed out of the instance ids.

    The method filter happens in SQL (idx_ddi_split_membership_method_split
    covers it) so rows for methods that never reach structural scoring are
    never pulled across into pandas at all."""
    placeholders = ",".join("?" for _ in methods)
    df = pd.read_sql_query(
        "SELECT DISTINCT ddi_id, instance_id_a, instance_id_b FROM ddi_split_membership "
        "WHERE instance_id_a IS NOT NULL AND instance_id_b IS NOT NULL "
        f"AND method IN ({placeholders})",
        conn,
        params=methods,
    )
    records = []
    for row in df.itertuples():
        pfam_id_a, uniprot_id_a, start_a, end_a = utils_struct.extract_data_from_instance_id(str(row.instance_id_a))
        pfam_id_b, uniprot_id_b, start_b, end_b = utils_struct.extract_data_from_instance_id(str(row.instance_id_b))
        records.append({
            "ddi_id": row.ddi_id,
            "instance_id_a": row.instance_id_a, "instance_id_b": row.instance_id_b,
            "pfam_id_a": pfam_id_a, "uniprot_id_a": uniprot_id_a, "start_a": start_a, "end_a": end_a,
            "pfam_id_b": pfam_id_b, "uniprot_id_b": uniprot_id_b, "start_b": start_b, "end_b": end_b,
        })
    columns = ["ddi_id", "instance_id_a", "instance_id_b",
               "pfam_id_a", "uniprot_id_a", "start_a", "end_a",
               "pfam_id_b", "uniprot_id_b", "start_b", "end_b"]
    return pd.DataFrame(records, columns=columns)


class BatchedInserter:
    """Buffers has_structure updates and flushes them with executemany() + a
    single commit per batch, instead of one UPDATE (and one commit per protein
    pair) per row. The UPDATE is an index lookup via
    idx_ddi_split_membership_triple (INIT_DOMAINSPLIT_DB)."""

    def __init__(self, conn: sqlite3.Connection, batch_size: int = INSERT_BATCH_SIZE):
        self.conn = conn
        self.batch_size = batch_size
        self._buf = []

    def add(self, ddi_id, instance_id_a, instance_id_b):
        self._buf.append((ddi_id, instance_id_a, instance_id_b))
        if len(self._buf) >= self.batch_size:
            self.flush()

    def flush(self):
        if not self._buf:
            return
        self.conn.executemany(
            "UPDATE ddi_split_membership SET has_structure = 1 "
            "WHERE ddi_id = ? AND instance_id_a = ? AND instance_id_b = ?",
            self._buf,
        )
        self.conn.commit()
        self._buf.clear()


def _join(values):
    return ",".join(sorted(set(values)))


def write_reports(db_path, methods, reasons, with_pae, with_iptm, coverage_path, to_model_dir):
    """structure_coverage.tsv and af3_to_model/{<method>,all_methods}.tsv from
    the output DB's membership rows of `methods`. reasons / with_pae /
    with_iptm are keyed by (instance_id_a, instance_id_b). Returns the
    coverage DataFrame."""
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    placeholders = ",".join("?" for _ in methods)
    df = pd.read_sql_query(
        "SELECT m.method, m.split, m.instance_id_a, m.instance_id_b, m.has_structure, d.negative "
        "FROM ddi_split_membership m JOIN domain_domain_interaction d ON d.id = m.ddi_id "
        "WHERE m.instance_id_a IS NOT NULL AND m.instance_id_b IS NOT NULL "
        f"AND m.method IN ({placeholders})", conn, params=methods)
    conn.close()
    # Same rule as write_run_manifest.py: the DDI's `negative` flag is the label.
    df["label"] = np.where(df["negative"].astype(bool), "negative", "positive")
    pairs = list(zip(df["instance_id_a"], df["instance_id_b"]))
    df["with_structure"] = df["has_structure"].astype(int)
    df["with_pae"] = [int(p in with_pae) for p in pairs]
    df["with_iptm"] = [int(p in with_iptm) for p in pairs]

    cov = (df.groupby(["method", "split", "label"])
             .agg(instance_pairs=("with_structure", "size"), with_structure=("with_structure", "sum"),
                  with_pae=("with_pae", "sum"), with_iptm=("with_iptm", "sum"))
             .reset_index())
    cov.insert(5, "fraction", (cov["with_structure"] / cov["instance_pairs"]).round(4))
    cov.to_csv(coverage_path, sep="\t", index=False)

    miss = df[df["with_structure"] == 0].copy()
    miss["reason"] = [reasons.get((a, b)) for a, b in zip(miss["instance_id_a"], miss["instance_id_b"])]
    if miss["reason"].isna().any():
        bad = miss[miss["reason"].isna()].iloc[0]
        raise RuntimeError(f"row without structure and without a recorded reason: "
                           f"{bad.method} {bad.split} {bad.instance_id_a}/{bad.instance_id_b}")
    ua = miss["instance_id_a"].map(lambda i: utils_struct.extract_data_from_instance_id(i)[1])
    ub = miss["instance_id_b"].map(lambda i: utils_struct.extract_data_from_instance_id(i)[1])
    miss["uniprot_id_a"], miss["uniprot_id_b"] = np.where(ua <= ub, ua, ub), np.where(ua <= ub, ub, ua)
    miss["instance_pair"] = miss["instance_id_a"] + "__" + miss["instance_id_b"]

    cols = ["uniprot_id_a", "uniprot_id_b", "reason", "n_instance_pairs", "labels", "splits"]
    os.makedirs(to_model_dir, exist_ok=True)

    def to_model(frame, extra):
        agg = dict(reason=("reason", _join), n_instance_pairs=("instance_pair", "nunique"),
                   labels=("label", _join), splits=("split", _join), **extra)
        out = (frame.groupby(["uniprot_id_a", "uniprot_id_b"]).agg(**agg).reset_index()
               if len(frame) else pd.DataFrame(columns=cols + list(extra)))
        return out[cols + list(extra)].sort_values(["uniprot_id_a", "uniprot_id_b"])

    for method in methods:
        to_model(miss[miss["method"] == method], {}).to_csv(
            os.path.join(to_model_dir, f"{method}.tsv"), sep="\t", index=False)
    to_model(miss, {"methods": ("method", _join)}).to_csv(
        os.path.join(to_model_dir, "all_methods.tsv"), sep="\t", index=False)
    return cov


def main():
    args = parse_args()
    n_workers = max(1, args.cpus if args.cpus else default_cpus())

    shutil.copy(args.db_in, args.db_out)
    conn = sqlite3.connect(args.db_out)
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA synchronous = NORMAL")
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]

    meta_index = build_meta_index(args.af_metadata)

    # Confirm metadata entries
    print(f"enrich_structural_af: built metadata index with {len(meta_index)} entries")

    con_ro = sqlite3.connect(f"file:{args.db_in}?mode=ro", uri=True)
    instances = fetch_relevant_instances(con_ro, methods)
    # Streamed and filtered: the protein table can hold every reviewed entry.
    needed = set(instances["uniprot_id_a"]) | set(instances["uniprot_id_b"])
    sequences = {u: seq for u, seq in con_ro.execute("SELECT uniprot_id, sequence FROM protein")
                 if u in needed}
    con_ro.close()

    # Unordered protein-pair key: (A, B) and (B, A) share one complex file, order lexicographically
    ua, ub = instances["uniprot_id_a"], instances["uniprot_id_b"]
    instances["pair_key"] = np.where(ua <= ub, ua + "|" + ub, ub + "|" + ua)

    n_pairs = instances["pair_key"].nunique()
    print(f"enrich_structural_af: {len(instances)} distinct instance-pairs across "
          f"{n_pairs} protein pairs to process (methods={methods}, "
          f"workers={n_workers})", flush=True)

    n_ok = n_missing_meta = n_missing_file = n_failed_slice = 0
    inserter = BatchedInserter(conn)
    reasons = {}                       # (instance_id_a, instance_id_b) -> why no structure
    with_pae, with_iptm = set(), set()  # (instance_id_a, instance_id_b) with that confidence
    conf_missing = []                  # confidence files that do not exist

    # Pass 1 (main process): pairs with a metadata entry become worker tasks;
    # the rest keep has_structure = 0.
    tasks = []
    for _, group in instances.groupby("pair_key", sort=False):
        first = group.iloc[0]
        entry = meta_index.get(frozenset((first.uniprot_id_a, first.uniprot_id_b)))
        if entry is None:
            n_missing_meta += len(group)
            for ia, ib in zip(group["instance_id_a"], group["instance_id_b"]):
                reasons[(ia, ib)] = "no_af_metadata_entry"
            continue
        path, file_uniprot_a, source = entry
        rows = list(group[ROW_COLS].itertuples(index=False, name=None))
        pair_seqs = {u: sequences.get(u) for u in (first.uniprot_id_a, first.uniprot_id_b)}
        tasks.append((path, source, file_uniprot_a, rows, pair_seqs, args.max_mismatch_rate))
    print(f"enrich_structural_af: {n_missing_meta} rows without metadata, "
          f"{len(tasks)} protein pairs to slice", flush=True)

    # Sort tasks by number of rows to slice, descending, so the longest-running
    # tasks are started first
    tasks.sort(key=lambda t: -len(t[3]))

    # Pass 2 (workers): load once, slice everything; main process writes results.
    t0 = time.time()
    n_done = 0
    with h5py.File(args.structures_h5, "w") as h5file:
        for res in run_tasks(tasks, n_workers):
            n_done += 1

            if res["unreadable_rows"]:
                print(f"  WARNING: {res['error']} -- {len(res['unreadable_rows'])} "
                      f"rows left without structure", flush=True)
                n_missing_file += len(res["unreadable_rows"])
                for _, ia, ib in res["unreadable_rows"]:
                    reasons[(ia, ib)] = res["unreadable_reason"]
            conf_missing += res["conf_missing"]

            for (group_key, dset_key), pdb_gz in res["slices"].items():
                grp = h5file.require_group(group_key)
                if dset_key not in grp:
                    # pdb_gz is already gzip-compressed, no h5 compression filter here
                    dset = grp.create_dataset(dset_key, data=np.frombuffer(pdb_gz, dtype="uint8"),
                                              compression=None)
                    attrs, pae, res_a, res_b = res["confidence"][(group_key, dset_key)]
                    dset.attrs.update(attrs)
                    pair = tuple(dset_key.split("__"))
                    if not np.isnan(attrs["iptm"]):
                        with_iptm.add(pair)
                    if pae is not None:
                        # Own top-level group: never collides with a {pfam_a}_{pfam_b} key.
                        pdset = h5file.require_group(f"pae/{group_key}").create_dataset(
                            dset_key, data=pae, compression="gzip")
                        pdset.attrs["residues_a"], pdset.attrs["residues_b"] = res_a, res_b
                        with_pae.add(pair)

            for ddi_id, ia, ib in res["ok_rows"]:
                inserter.add(ddi_id, ia, ib)
                n_ok += 1

            for ddi_id, ia, ib, err in res["failed"]:
                reasons[(ia, ib)] = "slice_failed"
                n_failed_slice += 1
                print(f"  WARNING: could not slice ddi_id={ddi_id} "
                      f"instances=({ia},{ib}): {err}", flush=True)

            if n_done % 500 == 0 or n_done == len(tasks):
                print(f"  progress: {n_done}/{len(tasks)} protein pairs "
                      f"({time.time() - t0:.0f}s)", flush=True)

    inserter.flush()
    # Hand off one self-contained file: a WAL-mode DB makes every later reader
    # (even mode=ro) create -wal/-shm sidecars next to it.
    conn.execute("PRAGMA journal_mode = DELETE")
    conn.close()
    print(f"enrich_structural_af: done -> sliced={n_ok} "
          f"missing_meta={n_missing_meta} missing_file={n_missing_file} "
          f"failed_slice={n_failed_slice}", flush=True)
    print(f"enrich_structural_af: confidence -> {len(with_pae)} sliced instance pairs with PAE, "
          f"{len(with_iptm)} with ipTM; {len(conf_missing)} confidence file(s) missing", flush=True)
    for missing in conf_missing[:20]:
        print(f"  WARNING: missing confidence file {missing}", flush=True)

    cov = write_reports(args.db_out, methods, reasons, with_pae, with_iptm,
                        args.coverage, args.af3_to_model_dir)
    for (method, label), g in cov.groupby(["method", "label"]):
        n_missing = int(g["instance_pairs"].sum() - g["with_structure"].sum())
        if n_missing:
            print(f"WARNING: enrich_structural_af: {n_missing}/{int(g['instance_pairs'].sum())} "
                  f"membership rows without structure (method={method}, label={label}); "
                  f"see af3_to_model/{method}.tsv", flush=True)

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")
        f.write(f"    biopython: {Bio.__version__}\n")
        f.write(f"    h5py: {h5py.__version__}\n")


if __name__ == "__main__":
    main()
