#!/usr/bin/env python3
"""Extract distinct Pfam accessions (without version suffix) from 3did SQLite DB.

Usage: 3did_extract_pfam_ids.py <3did.sqlite3> <pfam_ids.txt>
"""
import sqlite3
import sys

db_path, out_path = sys.argv[1], sys.argv[2]

con = sqlite3.connect(db_path)
cur = con.execute(
    "SELECT DISTINCT SUBSTR(Pfam_id, 1, INSTR(Pfam_id, '.') - 1) FROM Domain;"
)
with open(out_path, "w") as fh:
    for (pid,) in cur:
        if pid:
            fh.write(pid + "\n")
con.close()
