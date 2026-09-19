#!/usr/bin/env python3
"""Slice every AF3-predicted complex needed by ddi_split_membership into its
constituent domain-instance-pair structures, store them in structures.h5, and
record which (ddi_id, instance_id_a, instance_id_b) triples have a structure
in the `domain_structure` table.

Driven by ddi_split_membership, not by every domain x domain combination in
the database: only instance-pairs some (method, split) actually uses are
worth slicing. Grouped by protein pair so a complex predicted once is
resolved (CIF->PDB converted, if needed) once, even when it realizes several
DDIs.

Only the 'random' and 'minimal_leakage' methods get structural analysis
downstream, so fetch_relevant_instances restricts to those methods in SQL --
the other methods are what made the instance count explode (~1M -> ~250k),
and this way rows for methods that will never be scored are never even
pulled into pandas.

When af_metadata has no entry for a pair, or the referenced file is missing
(e.g. the AF3 run for it hasn't finished yet), the pair is marked is_mock=1
in `domain_structure` and no structure is written to structures.h5 -- 90%+
of pairs currently fall into this bucket, so mock rows are handled as a pure
metadata/DB write with no PDB parsing, slicing, or h5 I/O at all. Mock rows
are excluded from both build_scoring_matrix.py and score_ddi.py -- they
exist purely so the rest of the pipeline (H5 writing, DB inserts,
SUBSET_SPLIT_DB, etc.) can be exercised without waiting on every prediction
to land.
"""
import argparse as ap
import os
import shutil
import sqlite3
import sys
import tempfile

import Bio
import h5py
import numpy as np
import pandas as pd

import utils_struct

import glob
import zstandard

# Only these methods get structural analysis downstream (build_scoring_matrix.py,
# score_ddi.py) -- everything else in ddi_split_membership can be skipped before
# it ever leaves SQL.
STRUCTURAL_METHODS = ("random", "minimal_leakage")

# How many (ddi_id, instance_id_a, instance_id_b) rows to batch per
# executemany() call. Keeps a single Python list from growing unbounded on a
# very large run while still avoiding one round-trip per row.
INSERT_BATCH_SIZE = 5000


def parse_args():
    p = ap.ArgumentParser(description=__doc__, formatter_class=ap.RawDescriptionHelpFormatter)
    p.add_argument("--db_in", required=True)
    p.add_argument("--db_out", required=True)
    p.add_argument("--structures_h5", required=True)
    p.add_argument("--af_metadata", required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process_name", required=True)
    return p.parse_args()


def build_meta_index(meta_path: str) -> dict:
    """key: frozenset({uniprot_id_a, uniprot_id_b})
    value: (path, origin, file_uniprot_a) -- file_uniprot_a is whichever
    protein is chain A in the predicted file (af_metadata's own
    uniprot_id_a), so callers can tell whether their own (a, b) matches the
    file's orientation or is swapped."""
    meta = pd.read_csv(meta_path, sep=",", dtype=str)
    index = {}
    for _, row in meta.iterrows():
        key = frozenset({row.uniprot_id_a, row.uniprot_id_b})
        index[key] = (row.path, row.uniprot_id_a)  #  row.origin,
    return index


def resolve_structure(path: str, scratch_pdb: str) -> str | None:  # type: ignore
    """Return a plain .pdb path, decompressed from the AF3 model file under
    `path` (a directory containing '<model_entity_id>-model_v1.pdb.zst').
    Returns None if no matching file is found there."""
    candidates = glob.glob(os.path.join(path, "*-model_v1.pdb.zst"))
    if not candidates:
        return None
    if len(candidates) > 1:
        print(f"  WARNING: {len(candidates)} model files in {path}, using first: {candidates[0]}", flush=True)

    dctx = zstandard.ZstdDecompressor()
    with open(candidates[0], "rb") as fh_in, open(scratch_pdb, "wb") as fh_out:
        dctx.copy_stream(fh_in, fh_out)
    return scratch_pdb


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
    return pd.DataFrame(records)


def fetch_protein_sequences(conn) -> dict:
    return dict(conn.execute("SELECT uniprot_id, sequence FROM protein"))


class BatchedInserter:
    """Buffers domain_structure rows and flushes them with executemany()
    instead of one INSERT OR IGNORE per row. With ~90% of rows being mock
    (no PDB work at all), per-row round-trips to sqlite were themselves a
    meaningful share of the runtime."""

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
        self._buf.clear()


def main():
    args = parse_args()
    shutil.copy(args.db_in, args.db_out)
    conn = sqlite3.connect(args.db_out)
    # WAL + NORMAL synchronous: safe for a single-writer batch job like this
    # one (db_out is a fresh copy, nothing else touches it concurrently) and
    # noticeably faster than the default rollback-journal/FULL combo across
    # many small writes.
    conn.execute("PRAGMA journal_mode = WAL")
    conn.execute("PRAGMA synchronous = NORMAL")

    meta_index = build_meta_index(args.af_metadata)

    con_ro = sqlite3.connect(f"file:{args.db_in}?mode=ro", uri=True)
    instances = fetch_relevant_instances(con_ro)
    # protein_seqs = fetch_protein_sequences(con_ro)
    con_ro.close()

    n_pairs = instances[["uniprot_id_a", "uniprot_id_b"]].drop_duplicates().shape[0]
    print(f"enrich_structural_af: {len(instances)} distinct instance-pairs across "
          f"{n_pairs} protein pairs to process (methods={STRUCTURAL_METHODS})", flush=True)

    n_ok = n_mock = n_missing_meta = n_missing_file = n_missing_seq = 0
    inserter = BatchedInserter(conn)

    with tempfile.TemporaryDirectory() as scratch_dir, h5py.File(args.structures_h5, "w") as h5file:
        for (uniprot_a, uniprot_b), group in instances.groupby(["uniprot_id_a", "uniprot_id_b"]):
            is_mock = False
            entry = meta_index.get(frozenset((uniprot_a, uniprot_b)))
            scratch_pdb = os.path.join(scratch_dir, f"{uniprot_a}_{uniprot_b}.pdb")

            if entry is None:
                n_missing_meta += len(group)
                pdb_path = None
                is_mock, swapped = True, False
            else:
                path, file_uniprot_a = entry
                swapped = file_uniprot_a != uniprot_a
                pdb_path = resolve_structure(path, scratch_pdb)
                if pdb_path is None:
                    n_missing_file += len(group)
                    is_mock, swapped = True, False

            # Mock pairs skip PDB slicing entirely -- just record the rows.
            # This is the common case (~90% of pairs currently), so it's
            # worth keeping this branch a pure DB-buffer append with no
            # PDBParser/h5 work at all.
            if pdb_path is None:
                for row in group.itertuples():
                    inserter.add(row.ddi_id, row.instance_id_a, row.instance_id_b, is_mock)
                    n_mock += 1
                    n_ok += 1
                continue

            for row in group.itertuples():
                chain_a, chain_b = ("B", "A") if swapped else ("A", "B")
                try:
                    pdb_gz = utils_struct.ddi_pair_to_bytes(
                        pdb_path, chain_a, int(row.start_a), int(row.end_a),  # type: ignore
                        chain_b, int(row.start_b), int(row.end_b),  # type: ignore
                    )
                except Exception as exc:
                    print(f"  WARNING: could not slice ddi_id={row.ddi_id} "
                        f"instances=({row.instance_id_a},{row.instance_id_b}) "
                        f"from {pdb_path}: {exc}", flush=True)
                    continue

                group_key = f"{row.pfam_id_a}_{row.pfam_id_b}"
                dset_key = f"{row.instance_id_a}__{row.instance_id_b}"
                grp = h5file.require_group(group_key)
                if dset_key not in grp:
                    # pdb_gz is already gzip-compressed (see
                    # utils_struct.ddi_pair_to_bytes) -- no h5 compression
                    # filter here, or we'd be gzipping gzip for no benefit.
                    grp.create_dataset(dset_key, data=np.frombuffer(pdb_gz, dtype="uint8"),
                                        compression=None)

                inserter.add(row.ddi_id, row.instance_id_a, row.instance_id_b, is_mock=False)
                n_ok += 1

            inserter.flush()
            conn.commit()

    inserter.flush()
    conn.commit()
    conn.close()
    print(f"enrich_structural_af: done -> sliced={n_ok} (mock={n_mock}) "
          f"missing_meta={n_missing_meta} missing_file={n_missing_file} "
          f"missing_seq={n_missing_seq}", flush=True)

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")
        f.write(f"    biopython: {Bio.__version__}\n")
        f.write(f"    h5py: {h5py.__version__}\n")


if __name__ == "__main__":
    main()
