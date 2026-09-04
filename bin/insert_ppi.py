#!/usr/bin/env python3
"""Insert protein-protein interactions from STRING.

Reads a STRING links file (space-separated ``protein1 protein2 combined_score``)
and a gzipped UniProt ID mapping to translate STRING protein ids to UniProt
accessions, then inserts the scores into ``protein_protein_interaction``.

The mapping comes from PARSE_SWISSPROT's ``DR   STRING;`` lines rather than a
per-organism ``<ORG>_<taxid>_idmapping.dat.gz`` download; the format is the same
``uniprot \t id_type \t symbol`` triple. Nothing here is single-species: STRING
ids carry their own taxon prefix (``9606.ENSP…``) and are globally unique, so a
multi-species mapping and a multi-species links file join correctly with no
taxon bookkeeping at all.

**Streamed, not loaded.** The links file used to be one ``pd.read_csv`` -- fine
for the 9606 file, impossible for FETCH_STRING_LINKS' concatenation of every
taxon in the run. This reads it line by line and inserts in batches, so peak
memory is the mapping plus one batch regardless of how many species the run
covers.

Every dropped edge is counted, and counted **per taxon**: it used to be silent,
which meant a mapping that covered none of the links file -- the wrong organism,
or a taxon whose entries the ``--taxon-ids`` filter excluded -- produced zero PPI
rows and a green run. With many species in one file, a per-taxon breakdown is
what distinguishes "one species is misconfigured" from "the whole join is wrong".
"""

import argparse
import gzip
import shutil
import sqlite3
import sys
from collections import Counter, defaultdict
from typing import Dict

BATCH = 100_000

INSERT_SQL = (
    "INSERT OR IGNORE INTO protein_protein_interaction(protein_id_a, protein_id_b, score) "
    "VALUES (?, ?, ?)"
)


def load_uniprot_id_mapping(mapping_path: str, key_name: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with gzip.open(mapping_path, "rt") as f:
        for line in f:
            uniprot_id, id_type, symbol_name = line.strip().split("\t")
            if id_type.strip() == key_name:
                mapping[symbol_name.strip()] = uniprot_id.strip()
    return mapping


def open_links(path):
    """STRING links files are gzipped; a fixture may not be."""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def taxon_of(string_id: str) -> str:
    """``9606.ENSP00000269305`` -> ``9606``. STRING's own species-level taxid."""
    return string_id.split(".", 1)[0] if "." in string_id else "unknown"


def stream_edges(path, string_id_mapping, uniprot_to_pid, conn):
    """Insert every edge whose two endpoints resolve to a protein in the database.

    Returns ``(totals, per_taxon)``; ``per_taxon`` maps a STRING taxid to its own
    ``read``/``inserted``/``unmapped``/``not_in_db`` counter.
    """
    totals = Counter()
    per_taxon = defaultdict(Counter)
    batch = []

    with open_links(path) as fh:
        header = fh.readline()
        lines = fh if (header and header.startswith("protein1")) else _prepend(header, fh)
        for line in lines:
            parts = line.split()
            if len(parts) < 3:
                continue
            a, b, raw_score = parts[0], parts[1], parts[2]
            taxon = per_taxon[taxon_of(a)]
            totals["read"] += 1
            taxon["read"] += 1

            # Parse, don't pass through. This used to bind `parts[2]` straight
            # into the INSERT, and with a typeless `score` column SQLite stored
            # the string verbatim -- so every downstream consumer read TEXT '400'
            # and any numeric confidence filter (`score >= cutoff`) raised
            # instead of filtering. The column is REAL now, which would coerce a
            # bound string anyway; this is the second half of the same fix, and
            # it is where a genuinely non-numeric third column gets named.
            # Fatal, not skipped: STRING's links format is
            # `protein1 protein2 combined_score`, so a third field that is not a
            # number means this is not that file, and silently dropping the rows
            # would empty the interactome behind a green run.
            try:
                score = float(raw_score)
            except ValueError:
                raise SystemExit(
                    f"ERROR: {path}: edge {a} {b} (row {totals['read']} of the "
                    f"links file) has a non-numeric combined_score {raw_score!r}. "
                    "STRING links files are `protein1 protein2 combined_score` -- "
                    "check that --url_string / FETCH_STRING_LINKS' output is a "
                    "links file and not a different STRING download."
                ) from None

            uniprot_a = string_id_mapping.get(a)
            uniprot_b = string_id_mapping.get(b)
            if not uniprot_a or not uniprot_b:
                totals["unmapped"] += 1
                taxon["unmapped"] += 1
                continue
            pid_a = uniprot_to_pid.get(uniprot_a)
            pid_b = uniprot_to_pid.get(uniprot_b)
            if pid_a is None or pid_b is None:
                totals["not_in_db"] += 1
                taxon["not_in_db"] += 1
                continue

            batch.append((pid_a, pid_b, score))
            totals["inserted"] += 1
            taxon["inserted"] += 1
            if len(batch) >= BATCH:
                conn.executemany(INSERT_SQL, batch)
                batch.clear()

    if batch:
        conn.executemany(INSERT_SQL, batch)
    return totals, per_taxon


def _prepend(first, rest):
    """Put back a first line that turned out to be data, not a header."""
    if first:
        yield first
    yield from rest


def report(totals, per_taxon):
    """Print the funnel, then the per-taxon breakdown that localises a bad join.

    The two drop classes mean different things and only one of them is ordinary.
    ``not_in_db`` is expected and usually the large number: most STRING proteins
    carry no domain this run cares about. ``unmapped`` is the one to watch -- at
    scale it means the links file and the mapping disagree about which species
    they cover, and with several species in one file, which ones.
    """
    print(
        f"STRING edges: {totals['read']} read, {totals['inserted']} inserted, "
        f"{totals['unmapped']} with an unmapped STRING id, {totals['not_in_db']} whose "
        "proteins are not in this database",
        flush=True,
    )
    rows = sorted(per_taxon.items(), key=lambda kv: (-kv[1]["inserted"], kv[0]))
    shown = rows if len(rows) <= 25 else rows[:10]
    print(f"[ppi] per taxon ({len(rows)} in the links file{'' if len(rows) <= 25 else ', top 10 shown'}):", flush=True)
    for taxid, c in shown:
        print(
            f"[ppi]   {taxid}: read={c['read']} inserted={c['inserted']} "
            f"unmapped={c['unmapped']} not_in_db={c['not_in_db']}",
            flush=True,
        )
    dead = [t for t, c in rows if c["read"] and not c["inserted"]]
    if dead:
        print(
            f"[ppi] WARNING: {len(dead)} taxa contributed edges but none that resolved to "
            f"a protein in this database: {', '.join(dead[:20])}"
            + ("..." if len(dead) > 20 else ""),
            flush=True,
        )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-in", required=True)
    parser.add_argument("--string-ppi", required=True)
    parser.add_argument("--uniprot-id-mapping", required=True)
    parser.add_argument("--versions", required=True)
    parser.add_argument("--process-name", required=True)
    args = parser.parse_args()

    shutil.copy(args.db_in, "domainsplit.sqlite3")
    conn = sqlite3.connect("domainsplit.sqlite3")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    print("Inserting PPI", flush=True)

    uniprot_to_pid = dict(conn.execute("SELECT uniprot_id, id FROM protein").fetchall())
    print(f"Loaded {len(uniprot_to_pid)} protein ID mappings", flush=True)

    string_id_mapping = load_uniprot_id_mapping(args.uniprot_id_mapping, "STRING")
    print(f"Loaded {len(string_id_mapping)} STRING -> UniProt mappings", flush=True)

    totals, per_taxon = stream_edges(args.string_ppi, string_id_mapping, uniprot_to_pid, conn)
    conn.commit()
    conn.close()

    report(totals, per_taxon)

    if totals["read"] and not totals["inserted"]:
        raise SystemExit(
            f"ERROR: none of {totals['read']} STRING edges could be inserted "
            f"({totals['unmapped']} unmapped, {totals['not_in_db']} absent from the "
            "database). The STRING links file and the UniProt mapping most likely "
            "cover different species -- check --url_string, or FETCH_STRING_LINKS' "
            "string_taxa_report.tsv, against --instance_tier."
        )

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
