#!/usr/bin/env python3
"""Generate pooled domain embeddings for one model over one FASTA shard.

One model per invocation (`--model esm3|esmc|prott5|esm3_structure`). The
previous script loaded ESM3 *and* ESMC in a single task; splitting them gives
independent retries, a per-model batch size, and a shard's failure costs one
model instead of two.

There is only one output mode. Per-residue embeddings are gone: this pipeline
embeds the *cut domain sequence*, once per domain instance, because a DDI model
must not see protein context -- protein context correlates with the interaction
partner, so slicing a per-residue protein embedding smuggles protein identity
into a split that was partitioned on families.

H5 layout (the chunk contract, flat -- one model per file, so a `/model`
subgroup would be dead weight):

    <instance_id> -> (D,) float16

`instance_id` is whatever the FASTA header says, and ppi-splitting's
`sequences.fasta` is already keyed `{family}_{accession}_{start}_{end}`
(`fetch_domains.py`), so nothing here reconstructs a key.
EXPORT_DOMAIN_EMBEDDINGS re-keys these chunks to `{domain_id}/{instance_id}`.

Pooling is a masked mean, and the model families deliberately differ in what
the mask covers:

  - **ESM3 / ESMC / ESM3-structure**: mean over all `(BOS + seq + EOS)` tokens.
    Unchanged from the previous implementation, and deliberate -- ESM's
    terminal tokens carry sequence-level signal the model was trained to put
    there.
  - **ProtT5**: mean over the real residues only. T5 appends `</s>` and nothing
    else, and that position is a separator, not a summary; ProtBert/ProtT5's own
    embedding recipe drops it. Hence `attention_mask.sum(1) - 1`.

Structure embeddings (`--model esm3_structure`):
  - A separate model rather than a flag on `esm3`: structure conditioning only
    exists for ESM3 (ESMC / ProtT5 have no structure track), so keeping it a
    distinct `--model` value preserves "one model per task, one task per
    shard" instead of smuggling a second code path into the `esm3` branch.
  - Requires `--structures-mapping`, a CSV of `uniprot_id,pdb_id` (as written
    by `fetch_pdb_structures.py`). For each domain instance, the structure for
    its parent UniProt ID is fetched once via the ESM SDK
    (`ProteinChain.from_rcsb`) and cached, since many domain instances share a
    protein -- fetching per-domain would re-hit RCSB for the same chain
    repeatedly.
  - Alignment is best-effort: `fetch_pdb_structures.py` only supplies a
    `pdb_id`, not a UniProt-aligned residue map, so a domain's
    `start_pos`/`end_pos` (parsed from its instance id) are matched directly
    against the fetched chain's own `residue_index` numbers. A domain is only
    embedded with structure if every one of its residues resolves to a
    coordinate; anything else is skipped (absent from this model's H5, same
    as any other missing entry) and falls back to being covered by the
    plain `esm3` embedding instead.
"""

import argparse
import gc
import gzip
import os
import re
import sys

MODELS = ("esm3", "esmc", "prott5", "esm3_structure")

# ---------------------------------------------------------------------------
# Torch-free helpers. Deliberately importable without torch, h5py or any model
# weights, so tests/python/test_prott5_prep.py can assert the two things the
# ProtT5 branch is easy to get wrong without a GPU in the room.
# ---------------------------------------------------------------------------


def prepare_prott5_sequence(sequence: str) -> str:
    """ProtT5's expected input: rare residues mapped to X, residues whitespace-separated.

    `Rostlab/prot_t5_*` was trained on `[UZOB] -> X` substituted, space-separated
    residues; feeding it a bare sequence tokenizes to something else entirely.
    """
    return " ".join(re.sub(r"[UZOB]", "X", sequence.upper()))


def real_token_lengths(mask_sums) -> list:
    """Drop the trailing `</s>` from each ProtT5 attention-mask length.

    `attention_mask.sum(1)` counts real residues **plus** the appended `</s>`, so
    the pooled mean would otherwise average a separator token into the vector.
    Clamped at 1 so an empty sequence cannot produce a zero divisor.
    """
    return [max(int(total) - 1, 1) for total in mask_sums]


def masked_mean(emb, lengths):
    """Mean of `emb[i, :lengths[i], :]` over the sequence axis. `emb` is (B, L, D)."""
    import torch

    lengths_t = torch.tensor(lengths, device=emb.device).unsqueeze(1)   # (B,1)
    idx = torch.arange(emb.shape[1], device=emb.device).unsqueeze(0)    # (1,L)
    mask = (idx < lengths_t).to(emb.dtype).unsqueeze(-1)                # (B,L,1)
    summed = (emb * mask).sum(dim=1)                                    # (B,D)
    return summed / mask.sum(dim=1).clamp(min=1)                        # (B,D)


# ---------------------------------------------------------------------------
# HF cache / token
# ---------------------------------------------------------------------------


def _writable_dir(path: str) -> bool:
    """True if `path` exists (or can be created) and this user can write in it."""
    try:
        os.makedirs(path, exist_ok=True)
        probe = os.path.join(path, f".write_probe.{os.getpid()}")
        with open(probe, "w"):
            pass
        os.unlink(probe)
        return True
    except OSError:
        return False


def _setup_hf_cache(require_token: bool) -> None:
    # HF cache must be writable. A shared NFS dir (e.g. /nfs/scratch/hf_cache)
    # is preferred so parallel shards share one download, but it is often owned
    # by another user -- fall back to a task-local cache instead of dying.
    hf_home = os.environ.get("HF_HOME") or os.path.join(os.getcwd(), ".hf_cache")
    if not _writable_dir(hf_home):
        fallback = os.path.join(os.getcwd(), ".hf_cache")
        print(
            f"warning: HF_HOME={hf_home} is not writable; falling back to {fallback}",
            file=sys.stderr,
        )
        hf_home = fallback
        os.makedirs(hf_home, exist_ok=True)
    os.environ["HF_HOME"] = hf_home
    os.environ["HUGGINGFACE_HUB_CACHE"] = hf_home

    # ProtT5 is ungated, so a missing token is only fatal for ESM3 / ESMC.
    if require_token and not os.environ.get("HF_TOKEN"):
        raise RuntimeError(
            "HF_TOKEN not set. Required for gated ESM model download. "
            "Set via `nextflow secrets set HF_TOKEN <token>` and request "
            "access at https://huggingface.co/EvolutionaryScale/esm3-sm-open-v1"
        )
    # No huggingface_hub.login(): it persists the token to
    # $HF_HOME/stored_tokens, which fails on a shared cache dir and would leak
    # the token to everyone who can read it. huggingface_hub picks the token up
    # from the HF_TOKEN env var on its own for every authenticated request.


# ---------------------------------------------------------------------------
# FASTA
# ---------------------------------------------------------------------------


def _open_fasta(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _load_records(fasta_path: str, max_len: int):
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
    records.sort(key=lambda r: len(r[1]))
    return records


# ---------------------------------------------------------------------------
# Structure (esm3_structure only)
# ---------------------------------------------------------------------------


def _parse_domain_id(seq_id: str):
    """'{family}_{accession}_{start}_{end}' -> (uniprot_id, start, end) or None.

    Splits from the right so `family` may itself contain underscores.
    start/end are 1-based inclusive UniProt residue positions.
    """
    parts = seq_id.rsplit("_", 3)
    if len(parts) != 4:
        return None
    _family, uniprot_id, start_s, end_s = parts
    try:
        return uniprot_id, int(start_s), int(end_s)
    except ValueError:
        return None


def _load_structures_mapping(path: str) -> dict:
    """CSV of uniprot_id,pdb_id (as written by fetch_pdb_structures.py) -> dict."""
    import csv
    mapping = {}
    with open(path, "r") as f:
        for row in csv.DictReader(f):
            if row.get("pdb_id"):
                mapping[row["uniprot_id"]] = row["pdb_id"]
    return mapping


def _fetch_chain_coords(pdb_id: str):
    """Fetch one PDB chain and return {residue_number: (3,3 -> 37,3) coord row}.

    Keyed by the chain's own `residue_index` (not necessarily UniProt-numbered
    -- see module docstring), so lookups below are a direct dict access rather
    than an offset slice.
    """
    from esm.sdk.api import ESMProtein
    from esm.utils.structure.protein_chain import ProteinChain

    protein_chain = ProteinChain.from_rcsb(pdb_id.lower(), chain_id="detect")
    protein = ESMProtein.from_protein_chain(protein_chain)
    residue_numbers = protein_chain.residue_index.tolist()
    return {res_num: protein.coordinates[i] for i, res_num in enumerate(residue_numbers)}


class _StructureCache:
    """Per-uniprot_id chain cache so domains sharing a protein hit RCSB once."""

    def __init__(self, structure_mapping: dict):
        self._pdb_by_uniprot = structure_mapping
        self._chains = {}  # uniprot_id -> {residue_number: coord row} | None (failed)

    def coords_for_domain(self, uniprot_id: str, start_pos: int, end_pos: int):
        """Stacked (end-start+1, 37, 3) coords for a domain, or None if any
        residue in [start_pos, end_pos] is unresolved."""
        import numpy as np

        if uniprot_id not in self._chains:
            pdb_id = self._pdb_by_uniprot.get(uniprot_id)
            if not pdb_id:
                self._chains[uniprot_id] = None
            else:
                try:
                    self._chains[uniprot_id] = _fetch_chain_coords(pdb_id)
                except Exception as exc:
                    print(f"warn: could not fetch structure for {uniprot_id} (pdb {pdb_id}): {exc}", flush=True)
                    self._chains[uniprot_id] = None

        chain = self._chains[uniprot_id]
        if chain is None:
            return None
        rows = [chain.get(pos) for pos in range(start_pos, end_pos + 1)]
        if any(row is None for row in rows):
            return None
        import torch
        return torch.from_numpy(np.stack(rows).astype(np.float32))


# ---------------------------------------------------------------------------
# ESM (esm3 / esmc / esm3_structure)
# ---------------------------------------------------------------------------


def _encode_one(client, sequence: str, coordinates=None):
    """Tokenize a single sequence -> ESMProteinTensor (CPU-cheap)."""
    from esm.sdk.api import ESMProtein
    return client.encode(ESMProtein(sequence=sequence, coordinates=coordinates))


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
    from esm.utils.constants import esm3 as C
    import torch
    if len(tensors) == 1:
        bt = _BatchedESMProteinTensor.from_protein_tensor(tensors[0])
        return _move_batch_to_device(bt, device)
    max_len = max(t.sequence.shape[0] for t in tensors)
    bt = _BatchedESMProteinTensor.from_protein_tensor(tensors[0])

    # Pad EVERY present per-position track to max_len, not just `sequence`.
    # forward() indexes sequence against structure/secondary_structure/sasa
    # token-for-token (e.g. `sequence_tokens == C.SEQUENCE_BOS_TOKEN` is used
    # to masked_fill the structure track), so if only `sequence` is padded
    # while `structure` etc. keep tensors[0]'s original (shorter) length, you
    # get "size of tensor a (max_len) must match size of tensor b (orig_len)".
    track_pad_tokens = {
        "sequence": getattr(C, "SEQUENCE_PAD_TOKEN", 0),
        "structure": getattr(C, "STRUCTURE_PAD_TOKEN", 0),
        "secondary_structure": getattr(C, "SS8_PAD_TOKEN", 0),
        "sasa": getattr(C, "SASA_PAD_TOKEN", 0),
    }
    for name, pad_id in track_pad_tokens.items():
        if getattr(tensors[0], name, None) is None:
            continue
        padded = torch.full(
            (len(tensors), max_len), pad_id,
            dtype=tensors[0].sequence.dtype, device=device,
        )
        for i, t in enumerate(tensors):
            v = getattr(t, name, None)
            if v is None:
                continue
            src = v.to(device, non_blocking=True)
            padded[i, : src.shape[0]] = src
        setattr(bt, name, padded)

    return _move_batch_to_device(bt, device)


class _NullCtx:
    def __enter__(self):
        return self

    def __exit__(self, *_):
        return False


_DEVICE_DIAG_DONE = False


def _device_diag(client, batched, device):
    global _DEVICE_DIAG_DONE
    if _DEVICE_DIAG_DONE:
        return
    import torch
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
    seq_dev = getattr(getattr(batched, "sequence", None), "device", "?")
    print(
        f"DEVICE_DIAG: target={device} | batch={seq_dev} | "
        f"params={sorted(param_devs)} | buffers={sorted(buf_devs)} | first_embed={embed_dev}",
        flush=True,
    )
    _DEVICE_DIAG_DONE = True


def _write_pooled(out_h5, ids, pooled):
    import numpy as np
    arr = pooled.float().cpu().numpy().astype(np.float16)
    for i, sid in enumerate(ids):
        out_h5.create_dataset(sid, data=arr[i])


def _process_esm_batch(client, seqs, ids, out_h5, device, coords_list=None):
    """Encode + forward one ESM batch and write the pooled vectors.

    `coords_list`, if given, is a per-record (L,37,3) array or None, parallel
    to `seqs`/`ids` -- used only by `--model esm3_structure`.
    """
    import torch
    from esm.sdk.api import LogitsConfig

    if coords_list is None:
        coords_list = [None] * len(seqs)
    tensors = [_encode_one(client, s, c) for s, c in zip(seqs, coords_list)]
    batched = _stack_batch(tensors, device)
    # Per-item tokenized length, BOS+seq+EOS -- see the docstring for why the
    # terminal tokens stay in the ESM mean.
    lengths = [int(t.sequence.shape[0]) for t in tensors]
    _device_diag(client, batched, device)

    cfg = LogitsConfig(sequence=True, return_embeddings=True)
    with torch.inference_mode():
        autocast = (
            torch.autocast(device_type="cuda", dtype=torch.bfloat16)
            if torch.cuda.is_available() else _NullCtx()
        )
        with autocast:
            out = client.logits(batched, cfg)
        emb = out.embeddings  # (B, L_max, D)
        pooled = masked_mean(emb, lengths)

    _write_pooled(out_h5, ids, pooled)

    del tensors, batched, emb, pooled
    if torch.cuda.is_available():
        torch.cuda.empty_cache()


# ---------------------------------------------------------------------------
# ProtT5
# ---------------------------------------------------------------------------


def _load_prott5(model_name: str, device):
    from transformers import T5EncoderModel, T5Tokenizer

    tokenizer = T5Tokenizer.from_pretrained(model_name, do_lower_case=False, legacy=True)
    model = T5EncoderModel.from_pretrained(model_name)
    # The `-half` checkpoint is fp16 on disk; fp16 matmuls are unimplemented on
    # CPU, so a CPU fallback run has to be promoted to fp32.
    model = model.float() if device.type == "cpu" else model.half()
    return tokenizer, model.to(device).eval()


def _process_prott5_batch(bundle, seqs, ids, out_h5, device):
    import torch

    tokenizer, model = bundle
    encoded = tokenizer.batch_encode_plus(
        [prepare_prott5_sequence(s) for s in seqs],
        add_special_tokens=True,
        padding="longest",
        return_tensors="pt",
    )
    input_ids = encoded["input_ids"].to(device)
    attention_mask = encoded["attention_mask"].to(device)

    with torch.inference_mode():
        out = model(input_ids=input_ids, attention_mask=attention_mask)
        emb = out.last_hidden_state  # (B, L_max, D)
        lengths = real_token_lengths(attention_mask.sum(1).tolist())
        pooled = masked_mean(emb, lengths)

    _write_pooled(out_h5, ids, pooled)

    del encoded, input_ids, attention_mask, emb, pooled
    if torch.cuda.is_available():
        torch.cuda.empty_cache()


# ---------------------------------------------------------------------------
# Batch loop
# ---------------------------------------------------------------------------


def _run_batches(records, batch_size, process_batch):
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
                process_batch(seqs, ids)
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


def _write_versions(versions_path: str, process_name: str, model: str) -> None:
    import h5py
    import torch
    with open(versions_path, "w") as f:
        f.write(f'"{process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    torch: {torch.__version__}\n")
        f.write(f"    h5py: {h5py.__version__}\n")
        if model in ("esm3", "esmc", "esm3_structure"):
            import esm
            f.write(f"    esm: {esm.__version__}\n")
        else:
            import transformers
            f.write(f"    transformers: {transformers.__version__}\n")


def _resolve_device(require_gpu: bool):
    import torch

    if torch.cuda.is_available():
        return torch.device("cuda")
    if require_gpu:
        # Fail here, not later. A 1.4B model on CPU is not a slower route to the
        # same answer, it is a route to none: loading one per task on a shared
        # node OOM-killed five concurrent tasks (exit 137) ~90 s in, and the
        # retry did it again. The usual cause is the container not seeing the
        # GPU -- `--nv` missing for the *engine actually in use* (an apptainer
        # run picks up `apptainer.runOptions`, not `singularity.runOptions`), or
        # a SLURM request that binds the GPU to a job step rather than the batch
        # script. Both look identical from in here, so print what we can see.
        print("ERROR: --require-gpu was given but torch.cuda.is_available() is False.", file=sys.stderr)
        print(f"  CUDA_VISIBLE_DEVICES = {os.environ.get('CUDA_VISIBLE_DEVICES', '<unset>')!r}", file=sys.stderr)
        print(f"  SLURM_JOB_GPUS       = {os.environ.get('SLURM_JOB_GPUS', '<unset>')!r}", file=sys.stderr)
        print(f"  SLURM_STEP_GPUS      = {os.environ.get('SLURM_STEP_GPUS', '<unset>')!r}", file=sys.stderr)
        print(f"  torch {torch.__version__}, built for CUDA {torch.version.cuda}", file=sys.stderr)
        print(f"  /dev/nvidiactl present: {os.path.exists('/dev/nvidiactl')}", file=sys.stderr)
        print("  If /dev/nvidiactl is absent the container was started without --nv.", file=sys.stderr)
        return None
    print("warning: CUDA not available, falling back to CPU (very slow)", file=sys.stderr)
    return torch.device("cpu")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--output-h5", required=True)
    parser.add_argument("--versions", required=True)
    parser.add_argument("--process-name", required=True)
    parser.add_argument("--model", choices=MODELS, required=True)
    parser.add_argument("--prott5-model", default="Rostlab/prot_t5_xl_half_uniref50-enc",
                        help="HuggingFace id of the ProtT5 encoder (--model prott5 only)")
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--max-len", type=int, default=0, help="0 = no cap")
    parser.add_argument("--require-gpu", action="store_true",
                        help="abort instead of falling back to CPU when no GPU is visible")
    parser.add_argument(
        "--structures-mapping", default=None,
        help="CSV of uniprot_id,pdb_id (from fetch_pdb_structures.py). Required for --model esm3_structure.",
    )
    args = parser.parse_args()

    if args.model == "esm3_structure" and not args.structures_mapping:
        parser.error("--model esm3_structure requires --structures-mapping")

    _setup_hf_cache(require_token=args.model in ("esm3", "esmc", "esm3_structure"))

    import h5py
    import torch

    os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
    device = _resolve_device(args.require_gpu)
    if device is None:
        return 1

    records = _load_records(args.input_fasta, args.max_len)
    print(f"loaded {len(records)} records from {args.input_fasta}", flush=True)

    coords_by_id = None
    if args.model == "esm3_structure":
        structure_mapping = _load_structures_mapping(args.structures_mapping)
        cache = _StructureCache(structure_mapping)
        coords_by_id = {}
        kept = []
        for sid, seq in records:
            parsed = _parse_domain_id(sid)
            if parsed is None:
                print(f"skip {sid}: instance id doesn't parse as {{family}}_{{accession}}_{{start}}_{{end}}", flush=True)
                continue
            uniprot_id, start_pos, end_pos = parsed
            if end_pos - start_pos + 1 != len(seq):
                print(f"skip {sid}: domain span {start_pos}-{end_pos} doesn't match seq len {len(seq)}", flush=True)
                continue
            coords = cache.coords_for_domain(uniprot_id, start_pos, end_pos)
            if coords is None:
                continue
            coords_by_id[sid] = coords
            kept.append((sid, seq))
        print(f"{len(kept)}/{len(records)} domains have usable structure coordinates", flush=True)
        records = kept

    if not records:
        # Still write an empty H5 + versions so the export step doesn't crash.
        with h5py.File(args.output_h5, "w"):
            pass
        _write_versions(args.versions, args.process_name, args.model)
        return 0

    with h5py.File(args.output_h5, "w") as out_h5:
        print(f"loading {args.model}", flush=True)
        if args.model == "esm3":
            from esm.models.esm3 import ESM3
            client = ESM3.from_pretrained("esm3-open", device=device).to(device).eval()
            process = lambda seqs, ids: _process_esm_batch(client, seqs, ids, out_h5, device)  # noqa: E731
        elif args.model == "esm3_structure":
            from esm.models.esm3 import ESM3
            client = ESM3.from_pretrained("esm3-open", device=device).to(device).eval()
            process = lambda seqs, ids: _process_esm_batch(  # noqa: E731
                client, seqs, ids, out_h5, device,
                coords_list=[coords_by_id[i] for i in ids],
            )
        elif args.model == "esmc":
            from esm.models.esmc import ESMC
            client = ESMC.from_pretrained("esmc_600m", device=device).to(device).eval()
            process = lambda seqs, ids: _process_esm_batch(client, seqs, ids, out_h5, device)  # noqa: E731
        else:
            client = _load_prott5(args.prott5_model, device)
            process = lambda seqs, ids: _process_prott5_batch(client, seqs, ids, out_h5, device)  # noqa: E731

        _run_batches(records, args.batch_size, process)

        del client, process
        gc.collect()
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

    _write_versions(args.versions, args.process_name, args.model)
    return 0


if __name__ == "__main__":
    sys.exit(main())