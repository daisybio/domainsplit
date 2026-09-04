#!/usr/bin/env python3
"""Download STRING's per-organism links files for every species in the database,
keeping only the edges this run can use.

STRING publishes one links file per organism at a predictable URL::

    https://stringdb-downloads.org/download/protein.links.v<rel>/<taxid>.protein.links.v<rel>.txt.gz

The all-organism monolith (``protein.links.v12.0.txt.gz``) is **138 GB** and >99 %
of it is species this run has never heard of, so it is never used.

**The taxon list comes from the STRING ids, not from ``protein.taxon_id``.** A
STRING id is ``9606.ENSP00000269305`` and its prefix is STRING's own
species-level taxid, while UniProt's ``OX`` can be a strain-level id with no
STRING file behind it. Taking the prefix makes the join exact by construction.

Edges are filtered **here**, while streaming, rather than handed to
``insert_ppi.py`` in full: this step already holds the STRING -> UniProt map
restricted to the ``protein`` table, so an all-species concatenation would be
tens of GB of which a fraction of a percent survives the join. What
``insert_ppi.py`` receives is one small links file in STRING's own format.

A 404 for a taxon is counted and warned about, never fatal -- STRING does not
cover every NCBI taxid UniProt does. Running out of *every* taxon is fatal,
because a silently empty PPI table is the failure this whole path exists to
avoid.

Also reports the proteins that will never receive an interaction: UniProt does
not cross-reference STRING for unreviewed entries (measured on
``uniprot_trembl_human.dat.gz``: 4,616 entries carried **15** ``DR STRING;``
lines), so under a ``*_any_review_status`` tier the TrEMBL parents stay PPI-less
as a fact about UniProt, not as a pipeline bug.
"""

import argparse
import concurrent.futures
import gzip
import os
import sqlite3
import sys
import tempfile
from collections import Counter, defaultdict

import requests

URL_TEMPLATE = (
    "https://stringdb-downloads.org/download/protein.links.v{release}/"
    "{taxid}.protein.links.v{release}.txt.gz"
)

REPORT_COLUMNS = ["taxon_id", "proteins", "status", "edges_read", "edges_kept"]

CHUNK = 1 << 20


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--db", required=True, help="domainsplit SQLite; opened read-only")
    p.add_argument("--string-map", required=True, help="PARSE_SWISSPROT's uniprot_string_map.tsv.gz")
    p.add_argument("--out", required=True, help="filtered links file to write (.txt.gz)")
    p.add_argument("--report", required=True, help="per-taxon TSV to write")
    p.add_argument("--release", default="12.0", help="STRING release in the URL template")
    p.add_argument(
        "--min-proteins",
        type=int,
        default=1,
        help="skip a taxon contributing fewer than this many proteins, so one stray "
        "protein does not cost a 50 MB download",
    )
    p.add_argument("--cache-dir", default="", help="cache root; files land in <cache-dir>/string/")
    p.add_argument("--jobs", type=int, default=4, help="parallel downloads")
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def db_proteins(path):
    """``{uniprot_id: reviewed}`` for every row of ``protein``."""
    conn = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        rows = conn.execute("SELECT uniprot_id, reviewed FROM protein").fetchall()
    finally:
        conn.close()
    return {up: (rev or "") for up, rev in rows}


def load_string_map(path, proteins):
    """``({string_id: uniprot}, {taxid: {uniprot}})`` restricted to ``proteins``.

    The map file is ``uniprot \\t id_type \\t symbol``; only ``STRING`` rows count.
    """
    string_to_uniprot = {}
    per_taxon = defaultdict(set)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3 or parts[1].strip() != "STRING":
                continue
            uniprot, string_id = parts[0].strip(), parts[2].strip()
            if not string_id or uniprot not in proteins:
                continue
            string_to_uniprot[string_id] = uniprot
            per_taxon[string_id.split(".", 1)[0]].add(uniprot)
    return string_to_uniprot, per_taxon


def download(taxid, release, cache_dir, session):
    """``(path, status)``. ``path`` is None unless the file is on disk.

    ``status`` is ``cached``, ``downloaded``, ``404`` or ``error:<reason>``.
    """
    name = f"{taxid}.protein.links.v{release}.txt.gz"
    cached = os.path.join(cache_dir, name) if cache_dir else None
    if cached and os.path.exists(cached) and os.path.getsize(cached) > 0:
        return cached, "cached"

    url = URL_TEMPLATE.format(release=release, taxid=taxid)
    target_dir = cache_dir or os.getcwd()
    try:
        with session.get(url, stream=True, timeout=(30, 600)) as resp:
            if resp.status_code == 404:
                return None, "404"
            resp.raise_for_status()
            # Written to a temp file in the target directory and renamed, so a
            # killed task never leaves a truncated file behind for the next run
            # to read as a cache hit.
            fd, tmp = tempfile.mkstemp(dir=target_dir, prefix=f".{name}.", suffix=".part")
            try:
                with os.fdopen(fd, "wb") as out:
                    for chunk in resp.iter_content(CHUNK):
                        out.write(chunk)
                final = cached or os.path.join(target_dir, name)
                os.replace(tmp, final)
            except BaseException:
                if os.path.exists(tmp):
                    os.unlink(tmp)
                raise
    except requests.RequestException as exc:
        return None, f"error:{type(exc).__name__}"
    return final, "downloaded"


def stream_edges(path, string_to_uniprot, out):
    """Copy the edges of one links file whose both endpoints are wanted.

    STRING's format is space-separated ``protein1 protein2 combined_score`` with a
    header line. Kept verbatim, so ``insert_ppi.py`` reads the same shape it always
    did.
    """
    read = kept = 0
    with gzip.open(path, "rt") as fh:
        header = fh.readline()
        if header and not header.startswith("protein1"):
            # No header: the first line is data, so do not lose it.
            fh = _prepend(header, fh)
        for line in fh:
            parts = line.split()
            if len(parts) < 3:
                continue
            read += 1
            if parts[0] in string_to_uniprot and parts[1] in string_to_uniprot:
                out.write(f"{parts[0]} {parts[1]} {parts[2]}\n")
                kept += 1
    return read, kept


def _prepend(first, rest):
    yield first
    yield from rest


def main():
    args = parse_args()

    proteins = db_proteins(args.db)
    print(f"[string] {len(proteins)} proteins in the database", flush=True)

    string_to_uniprot, per_taxon = load_string_map(args.string_map, proteins)
    mapped = {up for ups in per_taxon.values() for up in ups}
    print(
        f"[string] {len(string_to_uniprot)} STRING ids map onto {len(mapped)} of them, "
        f"across {len(per_taxon)} taxa",
        flush=True,
    )

    # UniProt does not cross-reference STRING for unreviewed entries, so this
    # number is expected to be roughly the TrEMBL share of the universe. Reported
    # rather than investigated again in six months.
    unmapped = Counter(proteins[up] or "unknown" for up in proteins if up not in mapped)
    if unmapped:
        print(
            "[string] no STRING cross-reference, so no PPI edge, for "
            + ", ".join(f"{n} {flag}" for flag, n in sorted(unmapped.items()))
            + " proteins -- UniProt does not cross-reference STRING for unreviewed "
            "entries, and no links file changes that",
            flush=True,
        )

    if not per_taxon:
        raise SystemExit(
            "ERROR: none of this run's proteins carries a `DR STRING;` cross-reference, "
            "so there is no taxon to fetch links for and the PPI table would be empty. "
            "Check that PARSE_SWISSPROT's string map covers the same universe as the "
            "database, or pin one explicit links file with --url_string."
        )

    taxa = sorted(per_taxon, key=lambda t: (-len(per_taxon[t]), t))
    skipped = [t for t in taxa if len(per_taxon[t]) < args.min_proteins]
    taxa = [t for t in taxa if len(per_taxon[t]) >= args.min_proteins]
    if skipped:
        print(
            f"[string] {len(skipped)} taxa contribute fewer than {args.min_proteins} "
            "proteins and are skipped (--string_min_proteins_per_taxon)",
            flush=True,
        )
    if not taxa:
        raise SystemExit(
            f"ERROR: every one of the {len(per_taxon)} taxa contributes fewer than "
            f"{args.min_proteins} proteins, so nothing would be downloaded. Lower "
            "--string_min_proteins_per_taxon."
        )

    cache_dir = os.path.join(args.cache_dir, "string") if args.cache_dir else ""
    if cache_dir:
        os.makedirs(cache_dir, exist_ok=True)

    print(f"[string] fetching links for {len(taxa)} taxa (release {args.release})", flush=True)
    session = requests.Session()
    session.headers["User-Agent"] = "domainsplit-pipeline"
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(args.jobs, 1)) as pool:
        results = list(pool.map(lambda t: (t,) + download(t, args.release, cache_dir, session), taxa))

    rows, totals = [], Counter()
    with gzip.open(args.out, "wt", newline="") as out:
        out.write("protein1 protein2 combined_score\n")
        for taxid, path, status in results:
            read = kept = 0
            if path:
                read, kept = stream_edges(path, string_to_uniprot, out)
                if not cache_dir:
                    os.unlink(path)
            totals[status] += 1
            totals["edges_read"] += read
            totals["edges_kept"] += kept
            rows.append(
                {
                    "taxon_id": taxid,
                    "proteins": len(per_taxon[taxid]),
                    "status": status,
                    "edges_read": read,
                    "edges_kept": kept,
                }
            )

    for taxid in skipped:
        rows.append(
            {
                "taxon_id": taxid,
                "proteins": len(per_taxon[taxid]),
                "status": "skipped_min_proteins",
                "edges_read": 0,
                "edges_kept": 0,
            }
        )

    with open(args.report, "w", newline="") as fh:
        fh.write("\t".join(REPORT_COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in REPORT_COLUMNS) + "\n")

    missing = totals["404"] + sum(n for s, n in totals.items() if s.startswith("error:"))
    print(
        f"[string] taxa: {totals['downloaded']} downloaded, {totals['cached']} from cache, "
        f"{totals['404']} with no STRING file (404), "
        f"{sum(n for s, n in totals.items() if s.startswith('error:'))} failed; "
        f"edges: {totals['edges_read']} read, {totals['edges_kept']} kept",
        flush=True,
    )
    if missing:
        print(
            f"[string] WARNING: {missing} of {len(taxa)} taxa contributed no links file; "
            "their proteins get no interactions (see string_taxa_report.tsv)",
            flush=True,
        )
    if totals["downloaded"] + totals["cached"] == 0:
        raise SystemExit(
            f"ERROR: none of the {len(taxa)} taxa yielded a STRING links file "
            f"(release {args.release}). Either the release is wrong or the downloads "
            "are being blocked; see string_taxa_report.tsv for the per-taxon status."
        )
    if totals["edges_read"] and not totals["edges_kept"]:
        raise SystemExit(
            f"ERROR: {totals['edges_read']} STRING edges were read and none matched a "
            "protein in this database. The links files and the STRING id map disagree "
            "-- check the release against PARSE_SWISSPROT's `DR STRING;` ids."
        )

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    requests: {requests.__version__}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
