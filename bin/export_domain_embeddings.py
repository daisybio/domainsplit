#!/usr/bin/env python3
"""Re-key one model's embedding chunks into the layout the benchmark consumes.

Chunks arrive keyed flat by ``instance_id`` (see ``run_embeddings.py``). The
consumer -- ``daisybio-domainbenchmark``'s ``bin/features/embeddings.py`` -- reads

    h5[str(domain_id)][instance_key]

where ``domain_id`` is ``domain.id`` from the split database and ``instance_key``
is ``COALESCE(instance_id, 'r' || rowid)`` over ``domain_protein_map``. Only this
pipeline knows the mapping, so the re-key happens here, against the database, and
this step replaces the generic HDF5 join: one pass instead of a join followed by
a second full-file copy.

``domain.id`` is a **surrogate integer**. ``SUBSET_SPLIT_DB`` copies it verbatim
and ``PRUNE_UNREPRESENTED_DDIS`` deletes without renumbering, so one global file
is valid across every split database *of the same run* and silently wrong across
runs. The root attributes exist so a consumer can detect that rather than trust
it: ``domainsplit_run`` identifies the run that produced both artefacts.
"""

import argparse
import glob
import gzip
import sqlite3
import sys

import h5py
import numpy as np

ROOT_KEY_LAYOUT = "{domain_id}/{instance_id}"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--db", required=True, help="pruned domainsplit SQLite (read-only)")
    p.add_argument("--chunk-glob", default="chunk*",
                   help="glob for the staged per-shard HDF5 chunks")
    p.add_argument("--model", required=True, help="model name, recorded as a root attribute")
    p.add_argument("--output-h5", required=True)
    p.add_argument("--run-id", default="", help="value for the domainsplit_run root attribute")
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def _open_chunk(path):
    """Open a chunk, tolerating a gzipped one (the shards' own convention varies)."""
    try:
        opener = gzip.open(path, "rb")
        opener.peek(1)
        return h5py.File(opener, "r")
    except (OSError, gzip.BadGzipFile):
        return h5py.File(path, "r")


def read_instances(db_path):
    """``[(domain_id, instance_id, output_key)]`` for every domain instance.

    ``output_key`` is the benchmark's ``COALESCE(instance_id, 'r' || rowid)``.
    ``instance_id`` is what the chunks are keyed by, and is NULL only for a row no
    embedding can exist for -- the mapping is written once, by INGEST_INSTANCES,
    from the same ``instances.tsv`` the FASTA was keyed from.
    """
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    rows = conn.execute(
        "SELECT domain_id, instance_id, COALESCE(instance_id, 'r' || rowid) "
        "FROM domain_protein_map ORDER BY domain_id, rowid"
    ).fetchall()
    conn.close()
    return rows


def main() -> int:
    args = parse_args()

    chunks = sorted(glob.glob(args.chunk_glob))
    if not chunks:
        sys.exit(f"no chunks matched {args.chunk_glob!r}")
    print(f"[export] {len(chunks)} chunks for model {args.model}", flush=True)

    rows = read_instances(args.db)
    print(f"[export] {len(rows)} domain instances in {args.db}", flush=True)

    written = missing = no_instance_id = 0
    dims = set()
    domains = set()
    sample_missing = []

    with h5py.File(args.output_h5, "w") as out_h5:
        # Index the chunks once: which chunk holds which instance. The vectors are
        # (D,) fp16, so the index is tiny while the payload is not.
        located = {}
        for path in chunks:
            with _open_chunk(path) as chunk:
                for key in chunk.keys():
                    located[key] = path
        print(f"[export] {len(located)} embedded instances across the chunks", flush=True)

        by_chunk = {}
        for domain_id, instance_id, out_key in rows:
            if instance_id is None:
                no_instance_id += 1
                continue
            path = located.get(instance_id)
            if path is None:
                missing += 1
                if len(sample_missing) < 10:
                    sample_missing.append(instance_id)
                continue
            by_chunk.setdefault(path, []).append((domain_id, instance_id, out_key))

        for path in chunks:
            todo = by_chunk.get(path)
            if not todo:
                continue
            with _open_chunk(path) as chunk:
                for domain_id, instance_id, out_key in todo:
                    vector = np.asarray(chunk[instance_id], dtype=np.float16)
                    out_h5.create_dataset(f"{domain_id}/{out_key}", data=vector)
                    dims.add(int(vector.shape[-1]))
                    domains.add(domain_id)
                    written += 1

        out_h5.attrs["model"] = args.model
        out_h5.attrs["pooling"] = "mean"
        out_h5.attrs["dim"] = sorted(dims)[0] if len(dims) == 1 else -1
        out_h5.attrs["dtype"] = "float16"
        out_h5.attrs["key_layout"] = ROOT_KEY_LAYOUT
        out_h5.attrs["n_domains"] = len(domains)
        out_h5.attrs["n_instances"] = written
        out_h5.attrs["domainsplit_run"] = args.run_id

    print(f"[export] wrote {written} vectors under {len(domains)} domains -> {args.output_h5}", flush=True)
    if len(dims) > 1:
        print(f"[export] WARNING: mixed embedding dimensions {sorted(dims)}; dim attr set to -1", flush=True)
    if no_instance_id:
        print(f"[export] {no_instance_id} domain_protein_map rows have a NULL instance_id "
              "and cannot carry an embedding", flush=True)
    if missing:
        print(f"[export] WARNING: {missing} instances have no {args.model} embedding "
              f"(e.g. {', '.join(sample_missing)}) -- expected only for sequences dropped by "
              "--max-len or an OOM at batch size 1", flush=True)
    if written == 0:
        sys.exit("[export] no vector matched any domain instance -- the chunk keys and "
                 "domain_protein_map.instance_id have diverged")

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    h5py: {h5py.__version__}\n")
        f.write(f"    numpy: {np.__version__}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
