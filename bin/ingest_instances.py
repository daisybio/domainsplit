#!/usr/bin/env python3
"""Ingest ppi-splitting's domain instances into the domainsplit SQLite.

Reads the two files that describe the instance universe of the whole run --
``instances.tsv`` (one row per domain instance) and ``sequences.fasta`` (keyed by
instance id, holding the *domain* sequence, not the parent protein's) -- and
writes:

* ``domain``   -- one row per Pfam family seen,
* ``protein``  -- one row per parent UniProt accession (``sequence`` stays NULL;
  ENRICH_DDI_DATABASE fills it from the UniProt release),
* ``domain_protein_map`` -- one row per *instance*: ``domain_id, protein_id,
  domain_sequence, start_pos, end_pos, instance_id, clan, taxon_id``.

It also emits ``protein_domain_mapping.csv.gz`` in exactly the format
CREATE_PROTEIN_DOMAIN_MAPPING used to produce (``pfam_id,uniprot_id,start_pos,
end_pos,sequence``). Its one remaining consumer is
``insert_protein_sequences.py``, which reads only ``uniprot_id`` to pick the
UniProt records worth storing -- the ``sequence`` column is dead weight, kept
because the file is asserted line-for-line by
``tests/python/test_ingest_instances.py`` and an unused column is cheaper than a
churned interface. Domain embeddings are keyed by ``instance_id`` and read
ppi-splitting's ``sequences.fasta`` directly, so nothing reconstructs keys from
this file any more.

**Hard fail on any instance outside the configured taxon universe**
(``--taxon-ids``, from ``params.swissprot_taxon_ids``; empty accepts every
taxon).  An instance from outside it means the run's two taxon knobs disagree --
ppi-splitting sampled a species the protein universe was never parsed for -- and
carrying on would produce a database whose ``protein`` rows have no sequence and
no GO terms, silently.

The taxon is *kept*, not just checked: ``protein.taxon_id`` is what makes a
multi-species database stratifiable, and what tells a reader which proteins went
unenriched because their species had no STRING file.
"""

import argparse
import gzip
import sqlite3
import sys
from collections import Counter

from split_io import read_fasta, read_instances_tsv

MAPPING_HEADER = "pfam_id,uniprot_id,start_pos,end_pos,sequence\n"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True, help="domainsplit SQLite (modified in place)")
    p.add_argument("--instances", required=True, help="ppi-splitting instances.tsv")
    p.add_argument("--sequences", required=True, help="ppi-splitting sequences.fasta (keyed by instance id)")
    p.add_argument(
        "--taxon-ids",
        default="",
        help="comma-separated NCBI taxon ids the run's protein universe covers; "
        "an instance outside them is fatal. Empty accepts every taxon",
    )
    p.add_argument("--mapping-out", required=True, help="protein_domain_mapping.csv.gz to write")
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def load_instances(path, wanted_taxa):
    """Return ``(rows, stats)``; raises on any instance outside ``wanted_taxa``.

    ``wanted_taxa`` empty means every taxon is acceptable -- the any-species case.

    Duplicate ``instance_id``s are dropped -- ``instances.tsv`` is written once
    per run, but the DB's unique index on ``instance_id`` would fail loudly here
    rather than at the end of a long insert.
    """
    rows, seen, stats = [], set(), Counter()
    off_universe = []
    for row in read_instances_tsv(path):
        stats["read"] += 1
        taxon = (row["taxon_id"] or "").strip()
        if wanted_taxa and taxon not in wanted_taxa:
            off_universe.append((row["instance_id"], taxon or "<empty>"))
            continue
        iid = row["instance_id"].strip()
        if iid in seen:
            stats["duplicate_instance_id"] += 1
            continue
        seen.add(iid)
        rows.append(row)

    if off_universe:
        sample = ", ".join(f"{iid} (taxon {t})" for iid, t in off_universe[:10])
        raise SystemExit(
            f"ERROR: {len(off_universe)} of {stats['read']} instances are outside the "
            f"configured taxon universe ({', '.join(sorted(wanted_taxa))}): {sample}"
            + ("..." if len(off_universe) > 10 else "")
            + "\n--swissprot_taxon_ids decides which proteins exist and --instance_tier "
            "decides which of them may be sampled; these two disagree. Either widen "
            "--swissprot_taxon_ids or set --instance_tier human_only."
        )
    return rows, stats


def ingest(conn, rows, seqs):
    """Insert domains, proteins and instance-level mappings. Returns counters."""
    stats = Counter()

    families = sorted({r["family"].strip() for r in rows})
    conn.executemany("INSERT OR IGNORE INTO domain(pfam_id) VALUES (?)", [(f,) for f in families])
    stats["families"] = len(families)

    # One taxon per protein by construction -- an accession belongs to one entry --
    # so first-seen is the only value there is.
    protein_taxa = {}
    for r in rows:
        protein_taxa.setdefault(r["protein_id"].strip(), (r["taxon_id"] or "").strip())
    proteins = sorted(protein_taxa)
    conn.executemany(
        "INSERT OR IGNORE INTO protein(uniprot_id, taxon_id) VALUES (?, ?)",
        [(p, protein_taxa[p]) for p in proteins],
    )
    stats["proteins"] = len(proteins)
    stats["taxa"] = len(set(protein_taxa.values()))

    domain_id = {pfam: did for did, pfam in conn.execute("SELECT id, pfam_id FROM domain")}
    protein_id = {up: pid for pid, up in conn.execute("SELECT id, uniprot_id FROM protein")}

    mappings = []
    for row in rows:
        iid = row["instance_id"].strip()
        seq = seqs.get(iid)
        if seq is None:
            stats["no_sequence"] += 1
        mappings.append(
            (
                domain_id[row["family"].strip()],
                protein_id[row["protein_id"].strip()],
                seq,
                # int(), not the raw string: domain_protein_map.start_pos /
                # end_pos are INTEGER and part of the instance UNIQUE key, so a
                # TEXT binding would only match by affinity coercion. Binding the
                # integer keeps this writer honest independently of the schema.
                int(row["start"]),
                int(row["end"]),
                iid,
                (row["clan"] or "").strip() or None,
                row["taxon_id"].strip(),
            )
        )

    conn.executemany(
        "INSERT OR IGNORE INTO domain_protein_map"
        "(domain_id, protein_id, domain_sequence, start_pos, end_pos, instance_id, clan, taxon_id) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        mappings,
    )
    stats["instances"] = len(mappings)
    return stats


def write_mapping_csv(path, rows, seqs):
    """Write the CSV view INSERT_PROTEIN_SEQUENCES consumes, sorted for determinism."""
    written = 0
    with gzip.open(path, "wt", newline="") as out:
        out.write(MAPPING_HEADER)
        for row in sorted(rows, key=lambda r: (r["family"], r["protein_id"], int(r["start"]), int(r["end"]))):
            seq = seqs.get(row["instance_id"].strip())
            if not seq:
                continue
            out.write(f"{row['family']},{row['protein_id']},{row['start']},{row['end']},{seq}\n")
            written += 1
    return written


def main():
    args = parse_args()

    wanted_taxa = {t.strip() for t in args.taxon_ids.split(",") if t.strip()}
    print(
        f"[instances] taxon universe: {sorted(wanted_taxa) if wanted_taxa else 'none (all species)'}",
        flush=True,
    )
    rows, read_stats = load_instances(args.instances, wanted_taxa)
    seqs = read_fasta(args.sequences)
    print(f"[ingest] {read_stats['read']} instance rows, {len(rows)} kept, {len(seqs)} sequences", flush=True)
    if read_stats["duplicate_instance_id"]:
        print(f"[ingest] duplicate_instance_id = {read_stats['duplicate_instance_id']}", flush=True)

    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    stats = ingest(conn, rows, seqs)
    conn.commit()
    conn.close()

    written = write_mapping_csv(args.mapping_out, rows, seqs)

    for key in ("families", "proteins", "instances", "no_sequence"):
        print(f"[ingest] {key} = {stats[key]}", flush=True)
    print(f"[ingest] mapping_rows = {written}", flush=True)
    if stats["no_sequence"]:
        print(
            f"[ingest] WARNING: {stats['no_sequence']} instances have no sequence in "
            f"{args.sequences}; they are stored with a NULL domain_sequence and omitted "
            "from the mapping CSV, so they get no embedding.",
            flush=True,
        )

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
