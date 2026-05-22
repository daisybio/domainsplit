#!/usr/bin/env python3
"""Generate ESM3 + ESMC embeddings for a FASTA shard.

Two output modes:
  - per_residue : write (L+2, D) per sequence (BOS + seq + EOS, matches the
                  pre-batching pipeline output). Used for protein sequences
                  consumed downstream as per-residue tensors.
  - pooled      : mean across the (L+2) token dimension (BOS+EOS included, to
                  exactly match the legacy AVERAGE_POOL_EMBEDDINGS semantics).
                  Used for domain sequences.

Performance:
  - Batched logits() call via `_BatchedESMProteinTensor` (manual right-pad stack).
  - bf16 autocast around the model forward (autocast, NOT model.bfloat16(),
    since some ESM3 sub-modules are not bf16-safe).
  - Length bucketing: sort all records by length, batch contiguous groups.
  - On torch.OutOfMemoryError: halve the batch size and retry the same group.
    Cap halvings at 3; if the smallest batch still OOMs, fall back to
    per-record and skip records that still fail.
  - `--max-len` cap drops the long-tail quadratic-attention sequences entirely.
  - Storage dtype is float16 (downstream domainsplit only does np.array().dumps()).

H5 key contract (must match the legacy pipeline):
  per_residue: <seq_id>/esm3 -> (L+2, D), <seq_id>/esmc -> (L+2, D)
  pooled:      <seq_id>/esm3 -> (D,),     <seq_id>/esmc -> (D,)
"""

import argparse
import gc
import gzip
import os
import sys
from pathlib import Path


def _setup_hf_cache() -> None:
    # HF cache must be writable. If HF_HOME points at a shared NFS dir
    # (e.g. /nfs/scratch/hf_cache), use it; otherwise fall back to CWD.
    hf_home = os.environ.get("HF_HOME")
    if not hf_home:
        hf_home = os.path.join(os.getcwd(), ".hf_cache")
        os.environ["HF_HOME"] = hf_home
    os.environ.setdefault("HUGGINGFACE_HUB_CACHE", hf_home)
    try:
        os.makedirs(hf_home, exist_ok=True)
    except OSError as exc:
        print(f"warning: could not create HF_HOME={hf_home}: {exc}", file=sys.stderr)

    token = os.environ.get("HF_TOKEN")
    if not token:
        raise RuntimeError(
            "HF_TOKEN not set. Required for gated ESM model download. "
            "Set via `nextflow secrets set HF_TOKEN <token>` and request "
            "access at https://huggingface.co/EvolutionaryScale/esm3-sm-open-v1"
        )
    from huggingface_hub import login as _hf_login
    _hf_login(token=token, add_to_git_credential=False)


def _open_fasta(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _load_records(fasta_path: str, max_len: int, smoke_limit: int | None):
    from Bio import SeqIO

    with _open_fasta(fasta_path) as fh:
        records = []
        for rec in SeqIO.parse(fh, "fasta"):
            seq = str(rec.seq)
            if len(seq) == 0:
                continue
            if max_len > 0 and len(seq) > max_len:
                print(f"skip {rec.id}: len={len(seq)} > max_len={max_len}", flush=True)
                continue
            records.append((rec.id, seq))
            if smoke_limit and len(records) >= smoke_limit:
                break
    records.sort(key=lambda r: len(r[1]))
    return records


def _encode_one(client, sequence: str):
    """Tokenize a single sequence -> ESMProteinTensor (CPU-cheap)."""
    from esm.sdk.api import ESMProtein
    return client.encode(ESMProtein(sequence=sequence))


def _move_batch_to_device(bt, device):
    """Move every torch.Tensor field on `bt` to `device` (defensive against
    ESM SDK populating optional tracks like structure/sasa)."""
    import dataclasses
    import torch
    if dataclasses.is_dataclass(bt):
        attrs = [f.name for f in dataclasses.fields(bt)]
    else:
        attrs = [a for a in vars(bt).keys()]
    for name in attrs:
        v = getattr(bt, name, None)
        if isinstance(v, torch.Tensor) and v.device != device:
            setattr(bt, name, v.to(device, non_blocking=True))
    return bt


def _stack_batch(tensors, device):
    """Stack a list of ESMProteinTensor into a _BatchedESMProteinTensor on `device`."""
    import esm.sdk.api  # noqa: F401  prime to avoid circular import in esm 3.1.x
    from esm.utils.sampling import _BatchedESMProteinTensor
    import torch
    if len(tensors) == 1:
        bt = _BatchedESMProteinTensor.from_protein_tensor(tensors[0])
        return _move_batch_to_device(bt, device)
    max_len = max(t.sequence.shape[0] for t in tensors)
    pad_id = 0
    padded = torch.full((len(tensors), max_len), pad_id, dtype=tensors[0].sequence.dtype)
    for i, t in enumerate(tensors):
        padded[i, : t.sequence.shape[0]] = t.sequence
    bt = _BatchedESMProteinTensor.from_protein_tensor(tensors[0])
    bt.sequence = padded
    return _move_batch_to_device(bt, device)


def _seq_lengths(batched, tensors):
    """Per-item tokenized length (BOS+seq+EOS), as a list[int]."""
    return [int(t.sequence.shape[0]) for t in tensors]


def _process_batch(client, seqs, ids, mode: str, out_h5, model_key: str):
    """Encode + forward a batch, write outputs. Returns True if batch succeeded."""
    import numpy as np
    import torch
    from esm.sdk.api import LogitsConfig

    try:
        model_device = next(client.parameters()).device
    except StopIteration:
        model_device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")

    tensors = [_encode_one(client, s) for s in seqs]
    batched = _stack_batch(tensors, model_device)
    lengths = _seq_lengths(batched, tensors)

    cfg = LogitsConfig(sequence=True, return_embeddings=True)
    with torch.inference_mode():
        autocast = torch.autocast(device_type="cuda", dtype=torch.bfloat16) if torch.cuda.is_available() else _NullCtx()
        with autocast:
            out = client.logits(batched, cfg)
        emb = out.embeddings  # (B, L_max, D)

    if mode == "pooled":
        # Masked mean across the (BOS+seq+EOS) tokens, matching legacy AVERAGE_POOL_EMBEDDINGS.
        lengths_t = torch.tensor(lengths, device=emb.device).unsqueeze(1)  # (B,1)
        idx = torch.arange(emb.shape[1], device=emb.device).unsqueeze(0)   # (1,L)
        mask = (idx < lengths_t).to(emb.dtype).unsqueeze(-1)               # (B,L,1)
        summed = (emb * mask).sum(dim=1)                                   # (B,D)
        pooled = summed / mask.sum(dim=1).clamp(min=1)                     # (B,D)
        pooled = pooled.float().cpu().numpy().astype(np.float16)
        for i, sid in enumerate(ids):
            out_h5.create_dataset(f"{sid}/{model_key}", data=pooled[i])
    else:  # per_residue
        for i, sid in enumerate(ids):
            row = emb[i, : lengths[i], :].float().cpu().numpy().astype(np.float16)
            out_h5.create_dataset(f"{sid}/{model_key}", data=row)

    del tensors, batched, emb
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return True


class _NullCtx:
    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False


def _run_model(records, client, mode: str, batch_size: int, out_h5, model_key: str):
    """Iterate records in length-bucketed batches, with OOM halving."""
    import torch

    n = len(records)
    i = 0
    while i < n:
        # Try the current `batch_size`. On OOM, halve until 1, then skip.
        attempted = batch_size
        while attempted >= 1:
            end = min(i + attempted, n)
            ids = [r[0] for r in records[i:end]]
            seqs = [r[1] for r in records[i:end]]
            try:
                _process_batch(client, seqs, ids, mode, out_h5, model_key)
                i = end
                break
            except torch.OutOfMemoryError:
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
                print(
                    f"OOM at batch_size={attempted}, range [{i},{end}), "
                    f"len_range=[{len(seqs[0])},{len(seqs[-1])}]; halving",
                    flush=True,
                )
                attempted //= 2
        if attempted < 1:
            sid, seq = records[i]
            print(f"skip {sid}: cannot fit at batch_size=1 (len={len(seq)})", flush=True)
            i += 1


def _write_versions(versions_path: str, process_name: str) -> None:
    import esm
    import h5py
    import torch
    with open(versions_path, "w") as f:
        f.write(f'"{process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    esm: {esm.__version__}\n")
        f.write(f"    torch: {torch.__version__}\n")
        f.write(f"    h5py: {h5py.__version__}\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--output-h5", required=True)
    parser.add_argument("--versions", required=True)
    parser.add_argument("--process-name", required=True)
    parser.add_argument("--mode", choices=["per_residue", "pooled"], required=True)
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--max-len", type=int, default=0, help="0 = no cap")
    parser.add_argument("--smoke-limit", type=int, default=0, help="0 = no limit")
    args = parser.parse_args()

    _setup_hf_cache()

    import h5py
    import torch
    from esm.models.esm3 import ESM3
    from esm.models.esmc import ESMC

    os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
    device = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
    if not torch.cuda.is_available():
        print("warning: CUDA not available, falling back to CPU (very slow)", file=sys.stderr)

    smoke_limit = args.smoke_limit if args.smoke_limit > 0 else None
    records = _load_records(args.input_fasta, args.max_len, smoke_limit)
    print(f"loaded {len(records)} records from {args.input_fasta}", flush=True)
    if not records:
        # Still write an empty H5 + versions so downstream join doesn't crash.
        with h5py.File(args.output_h5, "w"):
            pass
        _write_versions(args.versions, args.process_name)
        return 0

    with h5py.File(args.output_h5, "w") as out_h5:
        print("loading ESM3", flush=True)
        client = ESM3.from_pretrained("esm3-open", device=device).to(device).eval()
        _run_model(records, client, args.mode, args.batch_size, out_h5, "esm3")
        del client
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        print("loading ESMC", flush=True)
        client = ESMC.from_pretrained("esmc_600m", device=device).to(device).eval()
        _run_model(records, client, args.mode, args.batch_size, out_h5, "esmc")
        del client
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    _write_versions(args.versions, args.process_name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
