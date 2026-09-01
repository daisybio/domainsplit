#!/usr/bin/env python3
"""Export the two files ppi-splitting needs from us, in its own formats.

Two modes, because the two exports sit on opposite sides of the Pfam fetch:

``--mode families``
    The *union* of every Pfam family the run will ever touch -- families in the
    DB (3did only at this point, by construction: external sources are parsed to
    TSVs and inserted last) plus families named in the external-source TSVs.
    One accession per line, which is what ``FETCH_DOMAIN_META`` takes.

    The union is fetched once so external families get domain instances too:
    ``BUILD_EXTERNAL_TEST`` needs concrete instance pairs for them, and
    ``domain_protein_map`` has to cover them for the enrichment steps.

``--mode ddis``
    The 3did positives as ppi-splitting's ``ppis`` CSV (``protein1,protein2``
    holding Pfam accessions in DDI mode), restricted to families that actually
    have an instance in the fetched ``instances.tsv``. That restriction is what
    keeps only 3did in the split population, and it drops the families
    ``--instance_tier human_reviewed`` stranded before BLAST ever sees them.
    ``PRUNE_UNREPRESENTED_DDIS`` removes the corresponding DDIs from the DB
    later; this only decides what gets split.

No negatives and no ``source`` column: negatives come back from ppi-splitting's
own samplers, and sources are re-joined from our DB at ingest.
"""

import argparse
import csv
import sqlite3
import sys

from ddi_db_utils import pfam_sort_key, split_sources
from external_ddi_tsv import read_external_tsv
from split_io import read_instances_tsv

PPIS_COLUMNS = ["protein1", "protein2"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--mode", required=True, choices=["families", "ddis"])
    p.add_argument("--db", required=True, help="domainsplit SQLite (read only)")
    p.add_argument(
        "--external-ddis",
        nargs="*",
        default=[],
        help="mode=families: normalized external-source TSVs whose families join the union",
    )
    p.add_argument(
        "--instances",
        help="mode=ddis: ppi-splitting's instances.tsv; families absent from it are dropped",
    )
    p.add_argument("--out", required=True, help="families.txt or ddis.csv")
    p.add_argument(
        "--source",
        default="3did",
        help="mode=ddis: only DDIs whose source list contains this source (default: 3did)",
    )
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    args = p.parse_args()
    if args.mode == "ddis" and not args.instances:
        p.error("--mode ddis requires --instances")
    return args


def db_families(conn):
    """Every Pfam accession referenced by a DDI in the database."""
    return {
        row[0]
        for row in conn.execute(
            "SELECT DISTINCT d.pfam_id FROM domain AS d "
            "WHERE d.id IN (SELECT domain_id_a FROM domain_domain_interaction "
            "               UNION SELECT domain_id_b FROM domain_domain_interaction)"
        )
    }


def source_ddis(conn, source):
    """``[(pfam_a, pfam_b)]`` for positive DDIs whose source list contains ``source``."""
    pairs = []
    for pfam_a, pfam_b, src in conn.execute(
        "SELECT da.pfam_id, db.pfam_id, ddi.source FROM domain_domain_interaction AS ddi "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b "
        "WHERE ddi.negative = 0"
    ):
        if source in split_sources(src):
            pairs.append((pfam_a, pfam_b))
    return pairs


def export_families(conn, external_paths, out_path):
    families = db_families(conn)
    n_db = len(families)
    for path in external_paths:
        for pfam_a, pfam_b, _negative, _source in read_external_tsv(path):
            families.add(pfam_a)
            families.add(pfam_b)
    ordered = sorted(families, key=pfam_sort_key)
    with open(out_path, "w") as fh:
        for pfam in ordered:
            fh.write(f"{pfam}\n")
    print(f"[export] families_in_db = {n_db}", flush=True)
    print(f"[export] families_external_only = {len(ordered) - n_db}", flush=True)
    print(f"[export] families_union = {len(ordered)}", flush=True)


def export_ddis(conn, source, instances_path, out_path):
    represented = {row["family"] for row in read_instances_tsv(instances_path)}
    pairs = source_ddis(conn, source)
    kept, dropped = [], 0
    for pfam_a, pfam_b in pairs:
        if pfam_a in represented and pfam_b in represented:
            kept.append((pfam_a, pfam_b))
        else:
            dropped += 1
    # Sorted so the rendered file -- and every task hash downstream of it -- does
    # not depend on SQLite's row order.
    kept.sort(key=lambda pair: (pfam_sort_key(pair[0]), pfam_sort_key(pair[1])))
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(PPIS_COLUMNS)
        writer.writerows(kept)
    print(f"[export] source = {source}", flush=True)
    print(f"[export] families_with_instances = {len(represented)}", flush=True)
    print(f"[export] ddis_total = {len(pairs)}", flush=True)
    print(f"[export] ddis_exported = {len(kept)}", flush=True)
    print(f"[export] ddis_dropped_no_instances = {dropped}", flush=True)


def main():
    args = parse_args()

    conn = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    if args.mode == "families":
        export_families(conn, args.external_ddis, args.out)
    else:
        export_ddis(conn, args.source, args.instances, args.out)
    conn.close()

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
