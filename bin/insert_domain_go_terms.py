#!/usr/bin/env python3
"""Insert domain GO term annotations from a pfam2go mapping file.

Parses the pfam2go text format to extract Pfam accession -> GO accession
mappings, then inserts them into the domain_go_terms table by joining on
the domain table's pfam_id.
"""

import argparse
import re
import shutil
import sqlite3
import sys

import pandas as pd


def load_pfam2go_dataframe(pfam2go_path: str) -> pd.DataFrame:
    with open(pfam2go_path) as fh:
        text = fh.read()
    matches = re.finditer(
        r"^Pfam:(?P<Pfam_accession>PF\d*)\s*(?P<Pfam_name>.*?)\s*>\s*GO:(?P<GO_name>.*?)\s*;\s*(?P<GO_accession>.*)\s*$",
        text,
        flags=re.MULTILINE,
    )
    return pd.DataFrame.from_records(match.groupdict() for match in matches)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-in", required=True)
    parser.add_argument("--pfam2go", required=True)
    parser.add_argument("--versions", required=True)
    parser.add_argument("--process-name", required=True)
    args = parser.parse_args()

    shutil.copy(args.db_in, "domainsplit.sqlite3")
    conn = sqlite3.connect("domainsplit.sqlite3")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    print("Inserting domain GO information", flush=True)
    pfam2go_df = load_pfam2go_dataframe(args.pfam2go)
    insert_iterator = (
        (row["GO_accession"], row["Pfam_accession"])
        for _, row in pfam2go_df.iterrows()
    )
    conn.executemany(
        """
        INSERT INTO domain_go_terms(domain_id, go_accession)
        SELECT domain.id as domain_id, ? AS go_accession
        FROM domain
        WHERE domain.pfam_id = ?;
        """,
        insert_iterator,
    )
    conn.commit()
    conn.close()

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    pandas: {pd.__version__}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
