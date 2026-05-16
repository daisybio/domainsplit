#!/usr/bin/env python3
"""Import a mysql2sqlite-converted SQL dump into a SQLite database.

Usage: 3did_import.py <dump.sql> <out.sqlite3>

mysql2sqlite prepends PRAGMA and BEGIN/END TRANSACTION lines; we handle
those ourselves so every row lands in a single bulk transaction instead of
one transaction per statement (which is what executescript() would give us,
since executescript auto-commits before each invocation).

Quote tracking: mysql2sqlite converts MySQL \\' -> '' and \\\\ -> \\ before
output, so the SQL contains no backslash escapes. We track '' (doubled
single-quote) and "" (doubled double-quote) -- both toggle the open-quote
flag twice in sequence, which correctly cancels out.
"""
import sqlite3
import sys

SKIP_PREFIXES = (
    "PRAGMA ",
    "BEGIN TRANSACTION",
    "END TRANSACTION",
    "COMMIT",
)

dump_path, db_path = sys.argv[1], sys.argv[2]

con = sqlite3.connect(db_path)
con.execute("PRAGMA synchronous = OFF")
con.execute("PRAGMA journal_mode = MEMORY")
con.execute("BEGIN")

buf = []
in_squote = False
in_dquote = False


def flush():
    stmt = "".join(buf).strip()
    if not stmt:
        return
    if any(stmt.upper().startswith(p) for p in SKIP_PREFIXES):
        return
    con.execute(stmt)


with open(dump_path, "r", encoding="utf-8", errors="replace") as fh:
    for line in fh:
        buf.append(line)
        for c in line:
            if c == "'" and not in_dquote:
                in_squote = not in_squote
            elif c == '"' and not in_squote:
                in_dquote = not in_dquote
        if not in_squote and not in_dquote and line.rstrip().endswith(";"):
            flush()
            buf = []

flush()
con.commit()
con.close()
