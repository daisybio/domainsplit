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
import re
import sqlite3
import sys

SKIP_PREFIXES = (
    "PRAGMA ",
    "BEGIN TRANSACTION",
    "END TRANSACTION",
    "COMMIT",
)

# insert_3did.py reads only these; the other ~17 tables of the dump are skipped.
# shortcut: statements split at line-final ';', so a multi-line trigger body still
# breaks the import (as it did before); 3did dumps carry no triggers today.
KEEP_TABLES = {"ddi1", "domain"}
TABLE_OF = re.compile(
    r"(?:INSERT(?:\s+OR\s+\w+)?\s+INTO|REPLACE\s+INTO|CREATE\s+TABLE(?:\s+IF\s+NOT\s+EXISTS)?"
    r"|DROP\s+TABLE(?:\s+IF\s+EXISTS)?|CREATE\s+(?:UNIQUE\s+)?INDEX\s+\S+\s+ON"
    r"|CREATE\s+TRIGGER\s+.*?\sON)\s+[`\"\[]?(\w+)",
    re.IGNORECASE | re.DOTALL,
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
    upper = stmt.upper()
    if any(upper.startswith(p) for p in SKIP_PREFIXES):
        return
    m = TABLE_OF.match(stmt)
    if m and m.group(1).lower() not in KEEP_TABLES:
        return
    # Some 3did dumps (notably the 2011 release) contain duplicate INSERT rows
    # that violate UNIQUE constraints once mysql2sqlite emits them as SQLite
    # tables. Rewrite plain INSERTs so the import skips duplicates instead of
    # aborting. Leave INSERT OR ... statements unchanged.
    if upper.startswith("INSERT INTO"):
        stmt = "INSERT OR IGNORE" + stmt[len("INSERT"):]
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
