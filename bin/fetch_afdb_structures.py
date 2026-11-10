#!/usr/bin/env python3
"""UniProt-numbered atom37 coordinates from AlphaFold DB, for `--model esm3_structure`.

One AFDB monomer model (`AF-<accession>-F1-model_v<version>.pdb`) per parent
protein of a domain instance. AFDB numbers residues by UniProt position, so a
residue's `residue_index` *is* its UniProt position and no sequence alignment is
needed -- but the model was predicted from whatever sequence version UniProt had
then, so every placed residue's one-letter code is still compared against this
run's UniProt sequence (the flat-file `SQ`, via PARSE_SWISSPROT's FASTA). A
protein above `--max-mismatch-rate` is rejected: its numbering no longer means
what the instance ids assume.

Residues with pLDDT below `--min-plddt` are written as NaN. ESM3 masks a
residue with non-finite backbone out of its geometric attention, and
run_embeddings.py also sets its structure token to <mask>, so a low-confidence
region conditions on sequence only instead of on a coordinate guess.

Output:
  --output-h5  one float32 dataset per accession, shape (L, 37, 3), L = the
               UniProt sequence length, NaN where unresolved. Root attrs record
               the AFDB version and pLDDT cutoff.
  --report     one row per accession: status, coverage, mismatch rate, detail.

Downloads land in `<cache-dir>/afdb-v<version>/` (gzip), so a second run is a
read. Runs on a CPU node; the GPU task only reads the H5.
"""

import argparse
import gzip
import io
import os
import sys
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor

URL_TEMPLATE = "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v{version}.pdb"
BACKBONE = (0, 1, 2)  # N, CA, C in atom37 order
REPORT_COLUMNS = ("uniprot_id", "status", "length", "placed", "mismatch_rate", "coverage", "detail")


# ---------------------------------------------------------------------------
# Pure helpers (numpy only), unit-tested in tests/python/test_afdb_structures.py
# ---------------------------------------------------------------------------


def place_on_uniprot(residue_index, chain_seq, atom37, plddt, uniprot_seq, min_plddt):
    """Place a chain's residues at their UniProt positions.

    `residue_index` holds 1-based UniProt positions (AFDB numbering). Returns
    `(coords, placed, mismatch_rate)`: `coords` is (len(uniprot_seq), 37, 3)
    float32, NaN where nothing was placed or pLDDT < `min_plddt`. A residue
    whose position falls outside the UniProt sequence counts as a mismatch,
    since that is the same symptom: the model's numbering is not this sequence's.
    """
    import numpy as np

    coords = np.full((len(uniprot_seq), 37, 3), np.nan, dtype=np.float32)
    placed = mismatches = 0
    for i, pos in enumerate(residue_index):
        placed += 1
        if not 1 <= pos <= len(uniprot_seq) or chain_seq[i] != uniprot_seq[pos - 1]:
            mismatches += 1
            continue
        if plddt[i] >= min_plddt:
            coords[pos - 1] = atom37[i]
    return coords, placed, (mismatches / placed if placed else 1.0)


def backbone_resolved(coords):
    """Boolean (L,): residues whose N, CA and C are all finite."""
    import numpy as np
    return np.isfinite(coords[:, BACKBONE, :]).all(axis=(1, 2))


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------


def accessions_from_instance_fasta(path):
    """Parent accessions of `{family}_{accession}_{start}_{end}` FASTA headers."""
    accs = set()
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for line in fh:
            if line.startswith(">"):
                parts = line[1:].split()[0].rsplit("_", 3)
                if len(parts) == 4:
                    accs.add(parts[1])
    return accs


def uniprot_sequences(path, wanted):
    """`>sp|ACC|ID ...` FASTA (PARSE_SWISSPROT) -> {accession: sequence} for `wanted`."""
    seqs, acc, buf = {}, None, []
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for line in fh:
            if line.startswith(">"):
                if acc in wanted:
                    seqs[acc] = "".join(buf)
                fields = line[1:].split("|")
                acc, buf = (fields[1] if len(fields) > 2 else line[1:].split()[0]), []
            else:
                buf.append(line.strip())
    if acc in wanted:
        seqs[acc] = "".join(buf)
    return seqs


def download(acc, version, cache_dir, session, retries=4, retry_wait=10.0):
    """`(path, status)`; status is cached / downloaded / 404 / error:<reason>.

    Same contract as fetch_string_links.download: a 404 is final at once, other
    failures are retried with doubling waits, and a file is only ever renamed
    into the cache complete.
    """
    import requests

    name = f"AF-{acc}-F1-model_v{version}.pdb.gz"
    final = os.path.join(cache_dir, name)
    if os.path.exists(final) and os.path.getsize(final) > 0:
        return final, "cached"
    status = None
    for attempt in range(max(retries, 1)):
        if attempt:
            time.sleep(retry_wait * 2 ** (attempt - 1))
        try:
            resp = session.get(URL_TEMPLATE.format(acc=acc, version=version), timeout=(30, 120))
            if resp.status_code == 404:
                return None, "404"
            resp.raise_for_status()
        except requests.RequestException as exc:
            status = f"error:{type(exc).__name__}"
            continue
        fd, tmp = tempfile.mkstemp(dir=cache_dir, prefix=f".{name}.", suffix=".part")
        try:
            with os.fdopen(fd, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb") as out:
                out.write(resp.content)
            os.replace(tmp, final)
        except BaseException:
            if os.path.exists(tmp):
                os.unlink(tmp)
            raise
        return final, "downloaded"
    return None, status


def parse_model(path):
    """AFDB PDB -> (residue_index, sequence, atom37, plddt) via ESM's own parser.

    ESM's parser, so the atom37 layout is by construction the one ESM3's
    structure encoder expects. `is_predicted` reads the B-factor as pLDDT.
    """
    from esm.utils.structure.protein_chain import ProteinChain

    with gzip.open(path, "rt") as fh:
        chain = ProteinChain.from_pdb(io.StringIO(fh.read()), chain_id="A", is_predicted=True)
    return chain.residue_index.tolist(), chain.sequence, chain.atom37_positions, chain.confidence


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--domain-fasta", required=True, help="instance-keyed domain FASTA (sequences.fasta)")
    p.add_argument("--uniprot-fasta", required=True, help="PARSE_SWISSPROT's uniprot_sequences.fasta.gz")
    p.add_argument("--output-h5", required=True)
    p.add_argument("--report", required=True)
    p.add_argument("--version", type=int, required=True, help="AFDB model version, e.g. 6")
    p.add_argument("--min-plddt", type=float, required=True)
    p.add_argument("--max-mismatch-rate", type=float, required=True)
    p.add_argument("--cache-dir", default="", help="cache root; files land in <cache-dir>/afdb-v<version>/")
    p.add_argument("--jobs", type=int, default=4)
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    args = p.parse_args()

    import h5py
    import numpy as np
    import requests

    accs = sorted(accessions_from_instance_fasta(args.domain_fasta))
    seqs = uniprot_sequences(args.uniprot_fasta, set(accs))
    print(f"[afdb] {len(accs)} parent proteins, {len(seqs)} with a UniProt sequence", flush=True)

    cache_dir = os.path.join(args.cache_dir or os.getcwd(), f"afdb-v{args.version}")
    os.makedirs(cache_dir, exist_ok=True)
    session = requests.Session()
    fetchable = [a for a in accs if a in seqs]
    with ThreadPoolExecutor(max_workers=max(args.jobs, 1)) as pool:
        fetched = dict(zip(fetchable, pool.map(
            lambda a: download(a, args.version, cache_dir, session), fetchable)))

    rows, n_written, n_errors = [], 0, 0
    with h5py.File(args.output_h5, "w") as h5:
        h5.attrs["afdb_version"] = args.version
        h5.attrs["min_plddt"] = args.min_plddt
        for acc in accs:
            row = dict.fromkeys(REPORT_COLUMNS, "")
            row["uniprot_id"] = acc
            rows.append(row)
            if acc not in seqs:
                row["status"] = "no_uniprot_sequence"
                continue
            row["length"] = len(seqs[acc])
            path, status = fetched[acc]
            if path is None:
                # AFDB has no F1 model for proteins over 2,700 residues, nor for
                # entries newer than its last sync: both are 404s.
                row["status"] = "not_in_afdb" if status == "404" else "download_failed"
                row["detail"] = status
                n_errors += status != "404"
                continue
            try:
                residue_index, chain_seq, atom37, plddt = parse_model(path)
            except Exception as exc:  # a truncated/corrupt cache file is per-protein, not fatal
                row["status"], row["detail"] = "parse_failed", f"{type(exc).__name__}: {exc}"[:200]
                continue
            coords, placed, rate = place_on_uniprot(
                residue_index, chain_seq, atom37, plddt, seqs[acc], args.min_plddt)
            row["placed"], row["mismatch_rate"] = placed, f"{rate:.4f}"
            if rate > args.max_mismatch_rate:
                row["status"] = "sequence_mismatch"
                row["detail"] = f"model has {placed} residues, UniProt {len(seqs[acc])}"
                continue
            row["coverage"] = f"{backbone_resolved(coords).mean():.4f}"
            row["status"] = "ok"
            h5.create_dataset(acc, data=coords, compression="gzip")
            n_written += 1

    with open(args.report, "w") as f:
        f.write("\t".join(REPORT_COLUMNS) + "\n")
        for row in rows:
            f.write("\t".join(str(row[c]) for c in REPORT_COLUMNS) + "\n")

    print(f"[afdb] {n_written}/{len(accs)} proteins written; {n_errors} download error(s)", flush=True)
    with open(args.versions, "w") as f:
        import esm
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    esm: {esm.__version__}\n")
        f.write(f"    h5py: {h5py.__version__}\n")
        f.write(f"    numpy: {np.__version__}\n")
        f.write(f"    afdb: v{args.version}\n")

    # Nothing usable at all is a broken network or a wrong --version, never data.
    if accs and n_written == 0:
        print("[afdb] ERROR: no protein got a structure; see the report", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
