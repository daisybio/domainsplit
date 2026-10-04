#!/usr/bin/env python3
"""Slice every AF3-predicted complex needed by ddi_split_membership into its
constituent domain-instance-pair structures, store them in structures.h5, and
record which (ddi_id, instance_id_a, instance_id_b) triples have a structure
via ddi_split_membership.is_mock.

Driven by ddi_split_membership, not by every domain x domain combination in
the database: only instance-pairs some (method, split) actually uses are
worth slicing. Grouped by UNORDERED protein pair so a complex predicted once
is loaded/parsed exactly once, even when it realizes several DDIs and even
when the pair shows up as (A, B) for some DDIs and (B, A) for others.

Only the 'random' and 'minimal_leakage' methods get structural analysis
downstream, so fetch_relevant_instances restricts to those methods in SQL --
the other methods are what made the instance count explode (~1M -> ~250k),
and this way rows for methods that will never be scored are never even
pulled into pandas.

When af_metadata has no entry for a pair, or the referenced file is missing
or unreadable (e.g. the AF3 run for it hasn't finished yet), the pair is
marked is_mock=1 and no structure is written to structures.h5 -- 90%+ of
pairs currently fall into this bucket, so mock rows are handled as a pure
metadata/DB write with no PDB parsing, slicing, or h5 I/O at all. Mock rows
are excluded from both build_scoring_matrix.py and score_ddi.py -- they
exist purely so the rest of the pipeline (H5 writing, DB inserts,
SUBSET_SPLIT_DB, etc.) can be exercised without waiting on every prediction
to land.

Parallelization: the expensive part (decompress/parse the complex once, then
slice + gzip every instance pair) runs in a pool of worker processes, one
task per protein pair. The main process is the only one that touches
structures.h5 and the sqlite DB (neither is safe for concurrent writers) and
consumes results as they arrive, with a bounded number of tasks in flight so
memory stays flat.
"""
import argparse as ap
import gzip
import io
import itertools
import multiprocessing
import os
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

import utils_struct

# Only these methods get structural analysis downstream (build_scoring_matrix.py,
# score_ddi.py) -- everything else in ddi_split_membership can be skipped before
# it ever leaves SQL.
STRUCTURAL_METHODS = ("random", "minimal_leakage")

# How many (ddi_id, instance_id_a, instance_id_b) rows to batch per
# executemany() + commit. Keeps a single Python list from growing unbounded on
# a very large run while still avoiding one round-trip/commit per row.
INSERT_BATCH_SIZE = 20000

# Temporary index so the per-row UPDATE ... WHERE ddi_id=? AND instance_id_a=?
# AND instance_id_b=? is an index lookup instead of a table scan. Dropped again
# at the end so the output DB schema is unchanged.
TMP_INDEX = "tmp_idx_dsm_enrich_triple"

# Columns handed to the workers, in this order.
ROW_COLS = ["ddi_id", "instance_id_a", "instance_id_b",
            "pfam_id_a", "uniprot_id_a", "start_a", "end_a",
            "pfam_id_b", "start_b", "end_b"]


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
    index = {}
    for row in meta.itertuples(index=False):
        key = frozenset({row.uniprot_id_a, row.uniprot_id_b})
        index[key] = (row.path, row.uniprot_id_a, row.source)
    return index


# ---------------------------------------------------------------------------
# Worker-side code (runs in child processes)
# ---------------------------------------------------------------------------

def load_structure(path: str, source: str):
    """Return a parsed Biopython Structure, fully in memory (no scratch files).
       source 'downloaded': `path` is a .pdb.zst AF3 model, decompressed on the fly
       source 'inferred':   `path` is a .cif file"""
    if source == "downloaded":
        dctx = zstandard.ZstdDecompressor()
        chunks = []
        with open(path, "rb") as fh:
            reader = dctx.stream_reader(fh)
            while True:
                chunk = reader.read(1 << 20)
                if not chunk:
                    break
                chunks.append(chunk)
        handle = io.StringIO(b"".join(chunks).decode("utf-8"))
        return PDBParser(QUIET=True).get_structure("complex", handle)
    if source == "inferred":
        return MMCIFParser(QUIET=True).get_structure("complex", path)
    raise ValueError(f"Unknown source {source} for path {path}")


def slice_pair(structure, chain_a: str, start_a: int, end_a: int,
               chain_b: str, start_b: int, end_b: int) -> bytes:
    """Same output as utils_struct.ddi_pair_to_bytes, but works on an
    already-parsed structure so the complex is parsed once per protein pair
    instead of once per instance pair. The domain from `chain_a` is written
    as chain A and the one from `chain_b` as chain B; the chain ids are
    restored afterwards so the structure can be reused for the next slice."""
    model = structure[0]
    ca, cb = model[chain_a], model[chain_b]
    # Two-step rename via temporary IDs: Biopython's Entity.id setter raises if
    # the new id already exists in the parent, so a direct A<->B swap would fail.
    ca.id, cb.id = "X", "Y"
    ca.id, cb.id = "A", "B"
    try:
        buf = io.StringIO()
        io_obj = PDBIO()
        io_obj.set_structure(structure)
        io_obj.save(buf, utils_struct.DdiPairSelect("A", start_a, end_a, "B", start_b, end_b))
        return gzip.compress(buf.getvalue().encode("utf-8"))
    finally:
        ca.id, cb.id = "X", "Y"
        ca.id, cb.id = chain_a, chain_b


def process_pair(task):
    """One protein pair: load the complex once, slice every instance pair.
    Never raises for data problems -- returns them so the main process can
    record them."""
    path, source, file_uniprot_a, rows = task
    try:
        structure = load_structure(path, source)
    except Exception as exc:  # missing / corrupt file -> treated as mock
        return {"unreadable_rows": [(r[0], r[1], r[2]) for r in rows],
                "error": f"{path}: {exc!r}",
                "slices": {}, "ok_rows": [], "failed": []}

    slices = {}    # (group_key, dset_key) -> gzip PDB bytes
    errors = {}    # (group_key, dset_key) -> error string
    ok_rows, failed = [], []
    for (ddi_id, ia, ib, pfam_a, uni_a, sa, ea, pfam_b, sb, eb) in rows:
        key = (f"{pfam_a}_{pfam_b}", f"{ia}__{ib}")
        if key not in slices and key not in errors:
            # Row's protein A is chain A in the file unless the file has the pair swapped.
            chain_a, chain_b = ("B", "A") if file_uniprot_a != uni_a else ("A", "B")
            try:
                slices[key] = slice_pair(structure, chain_a, int(sa), int(ea),
                                         chain_b, int(sb), int(eb))
            except Exception as exc:
                errors[key] = f"{exc!r} (from {path})"
        if key in errors:
            failed.append((ddi_id, ia, ib, errors[key]))
        else:
            ok_rows.append((ddi_id, ia, ib))
    return {"unreadable_rows": [], "error": None,
            "slices": slices, "ok_rows": ok_rows, "failed": failed}


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


def make_mock_structure(uniprot_a, uniprot_b, protein_seqs, scratch_dir):
    seq_a, seq_b = protein_seqs.get(uniprot_a), protein_seqs.get(uniprot_b)
    if not seq_a or not seq_b:
        print(f"WARNING: no sequence on file for ({uniprot_a}, {uniprot_b}) -- "
              f"cannot even mock a structure, skipping", flush=True)
        return None
    return utils_struct.mock_predict_complex(seq_a, seq_b, scratch_dir, f"{uniprot_a}_{uniprot_b}")


def fetch_relevant_instances(conn) -> pd.DataFrame:
    """Every (ddi_id, instance_id_a, instance_id_b) triple used by the
    structural methods in any split, deduplicated, with pfam/uniprot/
    coordinates parsed out of the instance ids.

    The method filter happens in SQL (idx_ddi_split_membership_method_split
    covers it) so rows for methods that never reach structural scoring are
    never pulled across into pandas at all -- this is what took the run from
    ~1M to ~250k candidate rows."""
    placeholders = ",".join("?" for _ in STRUCTURAL_METHODS)
    df = pd.read_sql_query(
        "SELECT DISTINCT ddi_id, instance_id_a, instance_id_b FROM ddi_split_membership "
        "WHERE instance_id_a IS NOT NULL AND instance_id_b IS NOT NULL "
        f"AND method IN ({placeholders})",
        conn,
        params=STRUCTURAL_METHODS,
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


def fetch_protein_sequences(conn) -> dict:
    return dict(conn.execute("SELECT uniprot_id, sequence FROM protein"))


class BatchedInserter:
    """Buffers is_mock updates and flushes them with executemany() + a single
    commit per batch, instead of one UPDATE (and one commit per protein pair)
    per row. With ~90% of rows being mock (no PDB work at all), per-row
    round-trips to sqlite were themselves a meaningful share of the runtime."""

    def __init__(self, conn: sqlite3.Connection, batch_size: int = INSERT_BATCH_SIZE):
        self.conn = conn
        self.batch_size = batch_size
        self._buf = []

    def add(self, ddi_id, instance_id_a, instance_id_b, is_mock: bool):
        self._buf.append((int(is_mock), ddi_id, instance_id_a, instance_id_b))
        if len(self._buf) >= self.batch_size:
            self.flush()

    def flush(self):
        if not self._buf:
            return
        self.conn.executemany(
            "UPDATE ddi_split_membership SET is_mock = ? "
            "WHERE ddi_id = ? AND instance_id_a = ? AND instance_id_b = ?",
            self._buf,
        )
        self.conn.commit()
        self._buf.clear()


def main():
    args = parse_args()
    n_workers = max(1, args.cpus if args.cpus else default_cpus())

    shutil.copy(args.db_in, args.db_out)
    conn = sqlite3.connect(args.db_out)
    # WAL + NORMAL synchronous: safe for a single-writer batch job like this
    # one (db_out is a fresh copy, nothing else touches it concurrently) and
    # noticeably faster than the default rollback-journal/FULL combo across
    # many small writes.
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA synchronous = NORMAL")
    conn.execute(f"CREATE INDEX IF NOT EXISTS {TMP_INDEX} "
                 "ON ddi_split_membership(ddi_id, instance_id_a, instance_id_b)")
    conn.commit()

    meta_index = build_meta_index(args.af_metadata)

    con_ro = sqlite3.connect(f"file:{args.db_in}?mode=ro", uri=True)
    instances = fetch_relevant_instances(con_ro)
    con_ro.close()

    # Unordered protein-pair key: (A, B) and (B, A) share one complex file, so
    # they must be one task (one load + parse), not two.
    ua, ub = instances["uniprot_id_a"], instances["uniprot_id_b"]
    instances["pair_key"] = np.where(ua <= ub, ua + "|" + ub, ub + "|" + ua)

    n_pairs = instances["pair_key"].nunique()
    print(f"enrich_structural_af: {len(instances)} distinct instance-pairs across "
          f"{n_pairs} protein pairs to process (methods={STRUCTURAL_METHODS}, "
          f"workers={n_workers})", flush=True)

    n_ok = n_mock = n_missing_meta = n_missing_file = n_failed_slice = 0
    inserter = BatchedInserter(conn)

    # Pass 1 (main process, no PDB work): mock pairs go straight to the DB
    # buffer; pairs with a metadata entry become worker tasks.
    tasks = []
    for _, group in instances.groupby("pair_key", sort=False):
        first = group.iloc[0]
        entry = meta_index.get(frozenset((first.uniprot_id_a, first.uniprot_id_b)))
        if entry is None:
            n_missing_meta += len(group)
            for row in group.itertuples():
                inserter.add(row.ddi_id, row.instance_id_a, row.instance_id_b, is_mock=True)
                n_mock += 1
            continue
        path, file_uniprot_a, source = entry
        rows = list(group[ROW_COLS].itertuples(index=False, name=None))
        tasks.append((path, source, file_uniprot_a, rows))
    inserter.flush()
    print(f"enrich_structural_af: {n_mock} rows mock (no metadata), "
          f"{len(tasks)} protein pairs to slice", flush=True)

    # Biggest tasks first so the last worker isn't stuck with a huge pair at the end.
    tasks.sort(key=lambda t: -len(t[3]))

    # Pass 2 (workers): load once, slice everything; main process writes results.
    t0 = time.time()
    n_done = 0
    with h5py.File(args.structures_h5, "w") as h5file:
        for res in run_tasks(tasks, n_workers):
            n_done += 1

            if res["unreadable_rows"]:
                print(f"  WARNING: {res['error']} -- marking {len(res['unreadable_rows'])} "
                      f"rows as mock", flush=True)
                for ddi_id, ia, ib in res["unreadable_rows"]:
                    inserter.add(ddi_id, ia, ib, is_mock=True)
                    n_mock += 1
                    n_missing_file += 1

            for (group_key, dset_key), pdb_gz in res["slices"].items():
                grp = h5file.require_group(group_key)
                if dset_key not in grp:
                    # pdb_gz is already gzip-compressed (see
                    # utils_struct.ddi_pair_to_bytes) -- no h5 compression
                    # filter here, or we'd be gzipping gzip for no benefit.
                    grp.create_dataset(dset_key, data=np.frombuffer(pdb_gz, dtype="uint8"),
                                       compression=None)

            for ddi_id, ia, ib in res["ok_rows"]:
                inserter.add(ddi_id, ia, ib, is_mock=False)
                n_ok += 1

            for ddi_id, ia, ib, err in res["failed"]:
                n_failed_slice += 1
                print(f"  WARNING: could not slice ddi_id={ddi_id} "
                      f"instances=({ia},{ib}): {err}", flush=True)

            if n_done % 500 == 0 or n_done == len(tasks):
                print(f"  progress: {n_done}/{len(tasks)} protein pairs "
                      f"({time.time() - t0:.0f}s)", flush=True)

    inserter.flush()
    conn.execute(f"DROP INDEX IF EXISTS {TMP_INDEX}")
    conn.commit()
    conn.close()
    print(f"enrich_structural_af: done -> sliced={n_ok} (mock={n_mock}) "
          f"missing_meta={n_missing_meta} missing_file={n_missing_file} "
          f"failed_slice={n_failed_slice}", flush=True)

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")
        f.write(f"    biopython: {Bio.__version__}\n")
        f.write(f"    h5py: {h5py.__version__}\n")


if __name__ == "__main__":
    main()