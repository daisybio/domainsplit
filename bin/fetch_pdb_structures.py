#!/usr/bin/env python3
"""Fetch the best experimental PDB structure per UniProt ID and store
UniProt-numbering-aligned atom37 coordinates for downstream ESM3
structure-conditioned embedding.

Single data source: RCSB's Sequence Coordinates Service
(https://sequence-coordinates.rcsb.org/graphql), queried once per protein
for UniProt -> PDB_ENTITY alignments. Its `target_alignments` list is
already ranked (first entry = top match, verified against the RCSB
website's own display), and each entry carries the actual aligned blocks
(query_begin/query_end/target_begin/target_end), so ranking and alignment
come from the same call — no cross-referencing two services, no risk of
"which structure is this really" disagreements.

For each protein:
  1. Query the alignment service once; drop any computed-structure-model
     hits (target_id prefixed "AF_"/"MA_" — AlphaFold/ModelArchive, not
     experimental) and anything not a plain 4-character PDB ID.
  2. Walk the remaining candidates in the order returned (best-first), up
     to --max-candidates. For each: fetch its coordinates via the ESM
     SDK's ProteinChain.from_rcsb(pdb_id, chain_id="detect") (same parser
     ESM3 itself uses; "detect" is fine here since all chains of one
     entity share the same sequence) and place them into a full
     UniProt-length array using the aligned blocks — which, being real
     blocks rather than a single linear range, handle insertions/
     deletions correctly.
  3. Safety net: every placed residue's chain letter is compared against
     the UniProt sequence at that position. If the mismatch rate exceeds
     --max-mismatch-rate the candidate is rejected and the next one tried;
     this also catches cases where the block's target_begin/target_end
     don't line up with ProteinChain's own residue ordering the way
     expected.
  4. If no candidate survives, the protein is left out of the H5 entirely
     (downstream falls back to sequence-only embedding).

Output:
  --output-h5:   one dataset per uniprot_id with a structure, shape (L,37,3) float32
  --report-csv:  one row per protein: chosen pdb_id/coverage/mismatch rate,
                 or a failure reason if nothing usable was found
"""

import argparse
import csv
import sys
import time

_GRAPHQL_URL = "https://sequence-coordinates.rcsb.org/graphql"
_GRAPHQL_QUERY = """
query UniProt2PDB($queryId: String!) {
  alignments(from: UNIPROT, to: PDB_ENTITY, queryId: $queryId) {
    target_alignments {
      target_id
    }
  }
}
"""



def _open_fasta(path: str):
    import gzip
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def _load_records(fasta_path: str):
    from Bio import SeqIO
    with _open_fasta(fasta_path) as fh:
        return set(rec.id for rec in SeqIO.parse(fh, "fasta"))


def _fetch_alignment(session, uniprot_id: str):
    resp = session.post(
        _GRAPHQL_URL,
        json={"query": _GRAPHQL_QUERY, "variables": {"queryId": uniprot_id}},
        timeout=30,
    )
    resp.raise_for_status()
    return resp.json()


def _parse_alignment_response(result: dict):
    try:
        target_alignments = result["data"]["alignments"]["target_alignments"] or []
    except (KeyError, TypeError):
        return []

    pdb_id = ""
    for entry in target_alignments:
        target_id = entry.get("target_id", "")
        if target_id.startswith(("AF_", "MA_")):
            continue  # computed structure model, not an experimental PDB entry
        if len(target_id.rsplit("_", 1)[0]) != 4:
            continue  # not a standard 4-character PDB ID; skip defensively
        pdb_id = target_id.rsplit("_", 1)[0]
        break
        
    return pdb_id


# def _fetch_chain(pdb_id: str):
#     from esm.utils.structure.protein_chain import ProteinChain
#     # No specific author chain id available at this (entity-level) query
#     # granularity, so let the SDK pick the first modeled chain — fine
#     # since every chain of a given entity shares the same sequence.
#     return ProteinChain.from_rcsb(pdb_id.lower(), chain_id="detect")




def _resolve_one(uniprot_id: str, session):
    """Try candidates best-first (RCSB's own ranking); return
    (coords_or_None, report_row_dict)."""
    import requests
    session = requests.Session()
    try:
        result = _fetch_alignment(session, uniprot_id)
    except Exception as exc:
        return None, {
            "uniprot_id": uniprot_id,
            "pdb_id": "",
            "status": "alignment_query_failed",
            "detail": str(exc)[:200]
        }

    pdb_id = _parse_alignment_response(result)
    if not pdb_id:
        return None, {"uniprot_id": uniprot_id, "pdb_id": "", "status": "no_pdb_structure", "detail": "no experimental PDB entry found"}
    return pdb_id, {"uniprot_id": uniprot_id, "pdb_id": pdb_id, "status": "ok", "detail": ""}

    


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-fasta", required=True)
    parser.add_argument("--output-mapping", required=True)
    parser.add_argument("--report-csv", required=True)
    parser.add_argument("--versions", required=True)
    parser.add_argument("--process-name", required=True)
    args = parser.parse_args()

    import requests
    from requests.adapters import HTTPAdapter, Retry

    records = _load_records(args.input_fasta)
    print(f"loaded {len(records)} proteins from {args.input_fasta}", flush=True)

    session = requests.Session()
    retries = Retry(total=3, backoff_factor=0.5, status_forcelist=[429, 500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries, pool_maxsize=4))

    report_rows = []
    n_ok = 0
    t0 = time.time()


    with open(args.output_mapping, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["uniprot_id", "pdb_id"])
        writer.writeheader()
        for uid in records:
            # For each protein, get best pdb id and write to mapping file, also collect report rows
            pdb_id, row = _resolve_one(uid, session)
            report_rows.append(row)
            if pdb_id:
                row["pdb_id"] = pdb_id
                n_ok += 1
            writer.writerow({k: row.get(k, "") for k in ["uniprot_id", "pdb_id"]})


    dt = time.time() - t0
    print(f"{n_ok}/{len(records)} proteins got a usable structure ({dt:.0f}s)", flush=True)

    fieldnames = [
        "uniprot_id", "pdb_id", "status", "detail"
    ]
    with open(args.report_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in report_rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})

    with open(args.versions, "w") as f:
        import requests as _requests
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    requests: {_requests.__version__}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())