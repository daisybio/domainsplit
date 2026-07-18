#!/usr/bin/env python3
"""Generate ESM3 + ESMC embeddings for a FASTA shard.

Two output modes:
  - per_residue : write (L+2, D) per sequence (BOS + seq + EOS). Used for
                  protein sequences consumed downstream as per-residue tensors.
  - pooled      : mean across all (L+2) tokens (BOS + seq + EOS included).
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


Structure embeddings:
  - If `--structures-h5` is provided, the ESM3 structure embeddings are also computed

H5 key contract (downstream pipeline expects this layout):
  per_residue: <seq_id>/esm3 -> (L+2, D), <seq_id>/esmc -> (L+2, D), <seq_id>/esm3_structure -> (L+2, D)
  pooled:      <seq_id>/esm3 -> (D,),     <seq_id>/esmc -> (D,),     <seq_id>/esm3_structure -> (D,)
"""

import argparse
import gc
import gzip
import os
import sys


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


# def _parse_domain_id(seq_id: str):
#     """'{pfam_id}_{uniprot_id}_{start}_{end}' -> (uniprot_id, start, end) or None.

#     Splits from the right so pfam_id may itself contain underscores.
#     start/end are 1-based inclusive UniProt residue positions.
#     """
#     parts = seq_id.rsplit("_", 3)
#     if len(parts) != 4:
#         return None
#     _pfam_id, uniprot_id, start_s, end_s = parts
#     try:
#         return uniprot_id, int(start_s), int(end_s)
#     except ValueError:
#         return None


# def _get_coords_for_record(seq_id: str, seq_len: int, mode: str, structures_h5):
#     """Look up (and validate) structure coordinates for one record.

#     Returns an (seq_len, 37, 3) float32 array, or None if unavailable /
#     invalid — callers fall back to sequence-only in that case.
#     """
#     import numpy as np

#     if mode == "pooled":
#         parsed = _parse_domain_id(seq_id)
#         if parsed is None:
#             return None
#         uniprot_id, start_pos, end_pos = parsed
#         if uniprot_id not in structures_h5:
#             return None
#         full = structures_h5[uniprot_id]
#         if start_pos < 1 or end_pos > full.shape[0] or start_pos > end_pos:
#             return None
#         coords = full[start_pos - 1: end_pos]
#     else:  # per_residue: seq_id is the uniprot_id itself
#         if seq_id not in structures_h5:
#             return None
#         coords = structures_h5[seq_id][:]

#     if coords.shape[0] != seq_len:
#         print(
#             f"warn: structure length mismatch for {seq_id} "
#             f"({coords.shape[0]} vs seq len {seq_len}); skipping entry",
#             flush=True,
#         )
#         return None
#     if np.isnan(coords).all():
#         return None
#     return coords.astype(np.float32)


def _fetch_structure(pdb_id: str):
    from esm.sdk.api import ESMProtein
    from esm.utils.structure.protein_chain import ProteinChain
    # Create a protein using a pdb format file from RCSB
    # Note: instead of the next two lines, we could use
    # protein_chain = ProteinChain.from_rcsb(pdb_id, chain_id)
    # but in future implementations, this function may use the mmcif file
    # which would throw off some indices later on in this notebook
    protein_chain = ProteinChain.from_rcsb(pdb_id.lower(), chain_id="detect")
    protein = ESMProtein.from_protein_chain(protein_chain)

    start = int(protein_chain.residue_index.min())
    end = int(protein_chain.residue_index.max())
    return protein.coordinates, start, end



def _attach_structures(records, mode: str, structure_mapping):
    """[(id, seq), ...] -> [(id, seq, coords_or_None), ...]."""
    if structure_mapping is None:
        return [(sid, seq, None, None) for sid, seq in records]
    out = []
    for sid, seq in records:
        # Extract uniprot_id from sid. For domain sequences, sid is of the form
        # "{pfam_id}_{uniprot_id}_{start}_{end}", so split from the right to allow underscores in pfam_id. For protein sequences, sid is the uniprot_id itself.
        if sid.count("_") >= 3:
            uniprot_id = sid.rsplit("_", 3)[1]
        else:
            uniprot_id = sid
            
        struct_seq = None
        coords = None

        if uniprot_id in structure_mapping:
            pdb_id = structure_mapping[uniprot_id]
            if not pdb_id or pdb_id == "":
                print(f"warn: no PDB ID for {sid} (uniprot {uniprot_id}); skipping entry", flush=True)
            else:
                try:
                    coords, start, end = _fetch_structure(pdb_id)
                    struct_seq = seq[start - 1 : end]
                    if coords.shape[0] != len(struct_seq):
                        print(
                            f"warn: structure length mismatch for {sid} "
                            f"(pdb {pdb_id} resolved range {start}-{end}: "
                            f"coords {coords.shape[0]} vs cropped seq {len(struct_seq)} "
                            f"[full seq len {len(seq)}]); skipping entry",
                            flush=True,

                        )
                        coords = None
                        struct_seq = None
                except Exception as exc:
                    print(f"warn: could not fetch structure for {sid} (pdb {pdb_id}): {exc}; skipping entry", flush=True)
        out.append((sid, seq, coords, struct_seq))
    return out



def _encode_one(client, sequence: str, coords = None):
    """Tokenize a single sequence -> ESMProteinTensor (CPU-cheap)."""
    from esm.sdk.api import ESMProtein

    # If seq = None, it was an empty structure sequence -> we want to skip it, so return an empty ESMProteinTensor
    if sequence is None:
        return client.encode(ESMProtein(sequence="")) # This should not make any problems downstream
    kwargs = {"sequence": sequence}
    if coords is not None:
        kwargs["coordinates"] = coords
    return client.encode(ESMProtein(**kwargs))


def _move_batch_to_device(bt, device):
    """Move every torch.Tensor field on `bt` to `device`.

    `_BatchedESMProteinTensor` is an attrs class (not a dataclass) and uses
    __slots__, so neither `dataclasses.fields` nor `vars()` enumerate its
    attributes. Prefer the SDK's `.to(device)` when present; fall back to a
    known-field list (defensive against the SDK populating optional tracks
    like structure/sasa)."""
    import torch
    to_fn = getattr(bt, "to", None)
    if callable(to_fn):
        try:
            result = to_fn(device)
            if result is not None:
                bt = result
        except Exception:
            pass
    known = (
        "sequence", "structure", "secondary_structure", "sasa",
        "function", "residue_annotations", "coordinates",
    )
    for name in known:
        v = getattr(bt, name, None)
        if isinstance(v, torch.Tensor) and v.device != device:
            try:
                setattr(bt, name, v.to(device, non_blocking=True))
            except (AttributeError, TypeError):
                pass
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
    # padded = torch.full(
    #     (len(tensors), max_len), pad_id,
    #     dtype=tensors[0].sequence.dtype, device=device,
    # )
    # for i, t in enumerate(tensors):
    #     src = t.sequence.to(device, non_blocking=True)
    #     padded[i, : src.shape[0]] = src
    # bt = _BatchedESMProteinTensor.from_protein_tensor(tensors[0])
    # bt.sequence = padded

    bt = _BatchedESMProteinTensor.from_protein_tensor(tensors[0])

    track_names = (
        "sequence", "structure", "secondary_structure", "sasa",
        "function", "residue_annotations", "coordinates",
    )


    for name in track_names:
        first = getattr(tensors[0], name, None)
        if not isinstance(first, torch.Tensor):
            continue  # track not populated for this batch (e.g. no coords given)
        padded = torch.full(
            (len(tensors), max_len) + tuple(first.shape[1:]),
            pad_id,
            dtype=first.dtype, device=device,
        )
        for i, t in enumerate(tensors):
            v = getattr(t, name, None)
            if not isinstance(v, torch.Tensor):
                # Shouldn't happen within one _run_model call (all records in
                # a call either all have coords or none do), but guard anyway.
                continue
            src = v.to(device, non_blocking=True)
            padded[i, : src.shape[0]] = src
        setattr(bt, name, padded)

    return _move_batch_to_device(bt, device)


def _seq_lengths(batched, tensors):
    """Per-item tokenized length (BOS+seq+EOS), as a list[int]."""
    return [int(t.sequence.shape[0]) for t in tensors]


_DEVICE_DIAG_DONE = False


def _process_batch(client, seqs, ids, coords, mode: str, out_h5, model_key: str, device):
    """Encode + forward a batch, write outputs. Returns True if batch succeeded."""
    import numpy as np
    import torch
    from esm.sdk.api import LogitsConfig

    tensors = [_encode_one(client, s, c) for s,c in zip(seqs, coords)]
    batched = _stack_batch(tensors, device)
    lengths = _seq_lengths(batched, tensors)

    global _DEVICE_DIAG_DONE
    if not _DEVICE_DIAG_DONE:
        try:
            param_devs = {str(p.device) for p in client.parameters()}
            buf_devs = {str(b.device) for b in client.buffers()}
        except Exception as e:
            param_devs = {f"err:{e}"}
            buf_devs = set()
        embed_dev = "?"
        for m in client.modules():
            if isinstance(m, torch.nn.Embedding):
                embed_dev = str(m.weight.device)
                break
        print(
            f"DEVICE_DIAG: target={device} | bt.sequence={batched.sequence.device} | "
            f"params={sorted(param_devs)} | buffers={sorted(buf_devs)} | first_embed={embed_dev}",
            flush=True,
        )
        _DEVICE_DIAG_DONE = True

    cfg = LogitsConfig(sequence=True, return_embeddings=True)
    with torch.inference_mode():
        autocast = torch.autocast(device_type="cuda", dtype=torch.bfloat16) if torch.cuda.is_available() else _NullCtx()
        with autocast:
            out = client.logits(batched, cfg)
        emb = out.embeddings  # (B, L_max, D)

    if mode == "pooled":
        # Masked mean across all (BOS+seq+EOS) tokens.
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


def _run_model(records, client, mode: str, batch_size: int, out_h5, model_key: str, device):
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
            coords = [r[2] for r in records[i:end]] if len(records[i]) > 2 else None
            try:
                _process_batch(client, seqs, ids, coords, mode, out_h5, model_key, device)
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
            sid, seq, _c = records[i]
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
    parser.add_argument(
        "--structures-mapping", default=None,
        help="Optional. Contains a CSV mapping of UniProt IDs to PDB structures. If provided, ESM3 structure embeddings will also be computed.",
    )
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

    structure_mapping = {}
    if args.structures_mapping:
        import csv
        with open(args.structures_mapping, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                uniprot_id = row["uniprot_id"]
                pdb_id = row["pdb_id"]
                structure_mapping[uniprot_id] = pdb_id

    records3 = _attach_structures(records, args.mode, structure_mapping)

    records_seq_only = [(sid, seq, None) for sid, seq, _c, _s in records3]
    records_with_struct = [(sid, struct_seq, coords) for sid, seq, coords, struct_seq in records3 if coords is not None]

    if records_with_struct:
        print(
            f"{len(records_with_struct)}/{len(records3)} records have usable "
            f"structure coordinates -> will also get an esm3_struct embedding",
            flush=True,
        )


    with h5py.File(args.output_h5, "w") as out_h5:
        print("loading ESM3", flush=True)
        client = ESM3.from_pretrained("esm3-open", device=device).to(device).eval()
        _run_model(records_seq_only, client, args.mode, args.batch_size, out_h5, "esm3", device)

        if records_with_struct:
            _run_model(records_with_struct, client, args.mode, args.batch_size, out_h5, "esm3_structure", device)

        
        del client
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        print("loading ESMC", flush=True)
        global _DEVICE_DIAG_DONE
        _DEVICE_DIAG_DONE = False
        client = ESMC.from_pretrained("esmc_600m", device=device).to(device).eval()
        _run_model(records_seq_only, client, args.mode, args.batch_size, out_h5, "esmc", device)
        del client
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()


    _write_versions(args.versions, args.process_name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
