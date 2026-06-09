#!/usr/bin/env python3
"""Infer positive DDIs from HIPPIE PPIs between two single-domain proteins.

A PPI contributes a positive DDI only when *both* interactors are reviewed human
proteins annotated with exactly one Pfam domain; the DDI is then the pair of
those two single domains.  Identifiers in the HIPPIE columns may be UniProt
accessions or entry names (e.g. ``AL1A1_HUMAN``) -- both are resolved via the
SwissProt map.  New domains are bulk-created so they get curated downstream.
"""

import argparse
import json
import sqlite3
import sys

from ddi_db_utils import count_source, ensure_domains, insert_ddis


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True)
    p.add_argument("--hippie", required=True)
    p.add_argument("--swissprot-map", required=True)
    p.add_argument("--min-score", type=float, required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def main():
    args = parse_args()

    with open(args.swissprot_map) as fh:
        smap = json.load(fh)
    accession_to_pfams = smap["accession_to_pfams"]
    name_to_accession = smap["name_to_accession"]

    # single-domain proteins: accession -> its one Pfam
    single_domain = {
        acc: pfams[0]
        for acc, pfams in accession_to_pfams.items()
        if len(pfams) == 1
    }
    print(f"[single_domain_ppi] single-domain proteins: {len(single_domain)}", flush=True)

    def resolve_pfam(token):
        """Return the single Pfam of ``token`` (accession or name), else None."""
        acc = token if token in accession_to_pfams else name_to_accession.get(token)
        if acc is None:
            return None
        return single_domain.get(acc)

    pairs = []
    n_rows = n_kept = n_unresolved = 0
    with open(args.hippie) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 5:
                continue
            n_rows += 1
            try:
                score = float(cols[4])
            except ValueError:
                continue
            if score < args.min_score:
                continue
            pfam_a = resolve_pfam(cols[0])
            pfam_b = resolve_pfam(cols[2])
            if pfam_a is None or pfam_b is None:
                n_unresolved += 1
                continue
            pairs.append((pfam_a, pfam_b))
            n_kept += 1

    print(f"[single_domain_ppi] hippie_rows={n_rows} score>= {args.min_score}: "
          f"single_domain_pairs={n_kept} unresolved_or_multi={n_unresolved}", flush=True)

    conn = sqlite3.connect(args.db)
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    pfams = {p for pair in pairs for p in pair}
    ensure_domains(conn, pfams)
    insert_ddis(conn, pairs, negative=False, source="single_domain_ppi")
    conn.commit()
    print(f"[single_domain_ppi] n_ddis_source = "
          f"{count_source(conn, 'single_domain_ppi')}", flush=True)
    conn.close()

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
