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
end_pos,sequence``), so ``generate_esm_embeddings`` and the ENRICH inserters keep
their current interface when CURATE_DOMAINS is deleted.

**Hard fail on any instance whose ``taxon_id`` is not 9606.**  ProtT5, STRING and
the UniProt idmapping steps downstream are all human-only; a non-human instance
here means ppi-splitting was not run with ``instance_tier = 'human_only'``, and
carrying on would produce a silently wrong database.
"""

import argparse
import gzip
import sqlite3
import sys
from collections import Counter

from split_io import read_fasta, read_instances_tsv

HUMAN_TAXON = "9606"

MAPPING_HEADER = "pfam_id,uniprot_id,start_pos,end_pos,sequence\n"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True, help="domainsplit SQLite (modified in place)")
    p.add_argument("--instances", required=True, help="ppi-splitting instances.tsv")
    p.add_argument("--sequences", required=True, help="ppi-splitting sequences.fasta (keyed by instance id)")
    p.add_argument("--mapping-out", required=True, help="protein_domain_mapping.csv.gz to write")
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def load_instances(path):
    """Return ``(rows, stats)``; raises on any non-human instance.

    Duplicate ``instance_id``s are dropped -- ``instances.tsv`` is written once
    per run, but the DB's unique index on ``instance_id`` would fail loudly here
    rather than at the end of a long insert.
    """
    rows, seen, stats = [], set(), Counter()
    non_human = []
    for row in read_instances_tsv(path):
        stats["read"] += 1
        taxon = (row["taxon_id"] or "").strip()
        if taxon != HUMAN_TAXON:
            non_human.append((row["instance_id"], taxon or "<empty>"))
            continue
        iid = row["instance_id"].strip()
        if iid in seen:
            stats["duplicate_instance_id"] += 1
            continue
        seen.add(iid)
        rows.append(row)

    if non_human:
        sample = ", ".join(f"{iid} (taxon {t})" for iid, t in non_human[:10])
        raise SystemExit(
            f"ERROR: {len(non_human)} of {stats['read']} instances are not human (taxon 9606): {sample}"
            + ("..." if len(non_human) > 10 else "")
            + "\nRun ppi-splitting with instance_tier = 'human_only'. Everything downstream of this "
            "step (ProtT5 embeddings, STRING PPI, UniProt idmapping) assumes human-only proteins."
        )
    return rows, stats


def ingest(conn, rows, seqs):
    """Insert domains, proteins and instance-level mappings. Returns counters."""
    stats = Counter()

    families = sorted({r["family"].strip() for r in rows})
    conn.executemany("INSERT OR IGNORE INTO domain(pfam_id) VALUES (?)", [(f,) for f in families])
    stats["families"] = len(families)

    proteins = sorted({r["protein_id"].strip() for r in rows})
    conn.executemany("INSERT OR IGNORE INTO protein(uniprot_id) VALUES (?)", [(p,) for p in proteins])
    stats["proteins"] = len(proteins)

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
                row["start"].strip(),
                row["end"].strip(),
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
    """Write the CSV view the ESM and ENRICH steps consume, sorted for determinism."""
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

    rows, read_stats = load_instances(args.instances)
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
