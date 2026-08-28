#!/usr/bin/env python3
"""Build the database for one ``(method, split)`` out of the enriched master.

``ddi_split_membership`` already records which DDIs belong to which split, so
this is a pure SQL filter: keep that split's DDIs, and keep exactly what a
surviving DDI reaches.  One module therefore replaces all three of the old
splitter modules -- the leakage-reduction logic lives in ppi-splitting, and what
is left here is bookkeeping.

**Copy-in, not copy-and-delete.**  This used to clone the whole master, ``DELETE``
almost all of it, then ``VACUUM`` -- which reads and writes the entire file twice
to produce a fraction of it, once per (method, split), 18 times. That was by far
the largest I/O in the pipeline when the master was ~99 % per-residue embedding
blobs. The blobs are gone -- embeddings are published as HDF5 and the master
carries none -- but the shape of the work has not changed: a clone-and-delete
still reads and writes every byte of the master 18 times to keep a fraction of
it, and protein sequences and the STRING network are not small. Now an empty
database is created, the master is ``ATTACH``ed **read-only**, and only the
surviving rows are inserted.

The schema is replayed from the master's own ``sqlite_master`` rather than
restated here, so there is no second source of truth to drift from
``INIT_DOMAINSPLIT_DB``; indexes are created after the inserts.

What "reaches" means, and it is the same reachability the delete cascade
implemented:

1. this split's DDIs;
2. the ``domain`` rows those DDIs reference -- and their ``domain_protein_map``
   and ``domain_go_terms`` rows;
3. the ``protein`` rows those mappings reference -- and their
   ``protein_go_terms`` rows, plus the ``protein_protein_interaction`` edges
   whose *both* endpoints survive (a cascade dropped an edge when either
   endpoint went, which is the same thing);
4. only this ``(method, split)``'s ``ddi_split_membership`` rows, so each
   published DB describes only itself.

Instances of a surviving family are kept whole, not narrowed to the pairs this
split happens to use: the split is a partition of *families*, so an instance of a
kept family belongs to this split by construction, and keeping them lets a
consumer form its own pairs.

Every table in the master must appear in ``FILTERS``. A new table added to
``INIT_DOMAINSPLIT_DB`` therefore fails here loudly, instead of being copied
whole into every split DB -- which is what the old delete-based version would
have done, silently publishing the rest of the run.
"""

import argparse
import os
import sqlite3
import sys
import urllib.parse

# Row-level reachability, per table. Each value is a WHERE clause evaluated
# against the attached master, referring to the keep_* temp tables built below.
FILTERS = {
    "domain_domain_interaction": "id IN (SELECT id FROM keep_ddi)",
    "domain": "id IN (SELECT id FROM keep_domain)",
    "domain_go_terms": "domain_id IN (SELECT id FROM keep_domain)",
    "domain_protein_map": "domain_id IN (SELECT id FROM keep_domain)",
    "protein": "id IN (SELECT id FROM keep_protein)",
    "protein_go_terms": "protein_id IN (SELECT id FROM keep_protein)",
    "protein_protein_interaction": (
        "protein_id_a IN (SELECT id FROM keep_protein) "
        "AND protein_id_b IN (SELECT id FROM keep_protein)"
    ),
    "ddi_split_membership": "method = :method AND split = :split",
}

# Reported before/after, in the order the old version printed them.
COUNTED = ["domain_domain_interaction", "domain", "protein", "domain_protein_map",
           "domain_go_terms", "protein_go_terms", "protein_protein_interaction"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--source", required=True, help="enriched master DB; opened read-only")
    p.add_argument("--out", required=True, help="split DB to create (must not exist)")
    p.add_argument("--method", required=True)
    p.add_argument("--split", required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def attach_source_readonly(conn, source):
    """``ATTACH`` the master as ``src``, read-only.

    Read-only is not a nicety: Nextflow stages the master as a symlink into the
    task directory, so a write would reach into another task's published output.
    The URI form is what carries ``mode=ro``, and it is honoured because the
    connection itself was opened with ``uri=True``.
    """
    uri = "file:" + urllib.parse.quote(os.path.abspath(source)) + "?mode=ro"
    conn.execute("ATTACH DATABASE ? AS src", (uri,))


def replay_schema(conn, kind):
    """Recreate the master's objects of one kind (``table`` / ``index`` / ...).

    Straight from ``src.sqlite_master``, so the split DB's schema is the master's
    by construction. ``sql IS NULL`` skips implicit indexes SQLite maintains for
    UNIQUE/PRIMARY KEY, which come back with their table.
    """
    rows = conn.execute(
        "SELECT name, sql FROM src.sqlite_master "
        "WHERE type = ? AND sql IS NOT NULL AND name NOT LIKE 'sqlite_%' "
        "ORDER BY rowid",
        (kind,),
    ).fetchall()
    for _name, sql in rows:
        conn.execute(sql)
    return len(rows)


def build_keep_tables(conn, method, split):
    """The three reachability sets, as temp tables keyed for fast lookup."""
    conn.execute("CREATE TEMP TABLE keep_ddi(id INTEGER PRIMARY KEY)")
    conn.execute(
        "INSERT INTO keep_ddi(id) SELECT DISTINCT ddi_id FROM src.ddi_split_membership "
        "WHERE method = :method AND split = :split AND ddi_id IS NOT NULL",
        {"method": method, "split": split},
    )

    conn.execute("CREATE TEMP TABLE keep_domain(id INTEGER PRIMARY KEY)")
    for column in ("domain_id_a", "domain_id_b"):
        conn.execute(
            f"INSERT OR IGNORE INTO keep_domain(id) SELECT DISTINCT {column} "
            "FROM src.domain_domain_interaction "
            f"WHERE id IN (SELECT id FROM keep_ddi) AND {column} IS NOT NULL"
        )

    conn.execute("CREATE TEMP TABLE keep_protein(id INTEGER PRIMARY KEY)")
    conn.execute(
        "INSERT OR IGNORE INTO keep_protein(id) SELECT DISTINCT protein_id "
        "FROM src.domain_protein_map "
        "WHERE domain_id IN (SELECT id FROM keep_domain) AND protein_id IS NOT NULL"
    )


def counts(conn, schema):
    return {t: conn.execute(f'SELECT COUNT(*) FROM {schema}."{t}"').fetchone()[0]
            for t in COUNTED}


def main():
    args = parse_args()

    if os.path.exists(args.out):
        sys.exit(f"[subset] {args.out} already exists; refusing to overwrite")

    # uri=True on the main connection is what lets ATTACH parse a file: URI.
    conn = sqlite3.connect(args.out, uri=True)
    # Left off for the build: the rows are copied table by table rather than in
    # dependency order, and the set is self-consistent by construction -- which
    # foreign_key_check verifies at the end rather than assuming.
    conn.execute("PRAGMA foreign_keys=OFF")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")
    attach_source_readonly(conn, args.source)

    tables = {r[0] for r in conn.execute(
        "SELECT name FROM src.sqlite_master WHERE type = 'table' AND name NOT LIKE 'sqlite_%'"
    )}
    missing = tables - set(FILTERS)
    if missing:
        sys.exit(f"[subset] no reachability rule for table(s) {sorted(missing)} -- add one to "
                 "FILTERS in bin/subset_split_db.py, or every split DB will silently omit them")

    before = counts(conn, "src")
    replay_schema(conn, "table")
    build_keep_tables(conn, args.method, args.split)

    kept = conn.execute("SELECT COUNT(*) FROM keep_ddi").fetchone()[0]
    if kept == 0:
        print(
            f"[subset] WARNING: no ddi_split_membership rows for {args.method}/{args.split}; "
            "the emitted database will be empty.",
            flush=True,
        )

    params = {"method": args.method, "split": args.split}
    for table in sorted(tables):
        conn.execute(
            f'INSERT INTO main."{table}" SELECT * FROM src."{table}" WHERE {FILTERS[table]}',
            params,
        )
    # After the bulk insert, not before: index maintenance per row would be paid
    # for nothing.
    n_index = replay_schema(conn, "index")
    for kind in ("view", "trigger"):
        replay_schema(conn, kind)
    conn.commit()

    after = counts(conn, "main")
    print(f"[subset] {args.method}/{args.split}: {kept} DDIs kept, {n_index} indexes rebuilt",
          flush=True)
    for table in COUNTED:
        print(f"[subset] {table} = {after[table]} (was {before[table]})", flush=True)

    conn.execute("PRAGMA foreign_keys=ON")
    violations = conn.execute("PRAGMA main.foreign_key_check").fetchall()
    if violations:
        sys.exit(f"[subset] {len(violations)} foreign key violation(s) in {args.out}: "
                 f"{violations[:5]}")
    conn.execute("DETACH DATABASE src")
    conn.close()

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
