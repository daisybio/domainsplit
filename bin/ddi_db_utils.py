#!/usr/bin/env python3
"""Shared helpers for inserting DDIs into the domainsplit SQLite.

Every ``INSERT_<source>`` step uses these so all sources are handled uniformly:
first bulk-create any missing ``domain`` rows for the Pfam IDs it references,
then insert its DDIs.  Pairs are order-normalised by Pfam accession number (see
:func:`pfam_sort_key`) so the stored ``(domain_id_a, domain_id_b)`` order is
stable -- a pair is deduplicated regardless of the order it is supplied in and
regardless of which source inserted it first or in what order domains were
created -- and ``INSERT OR IGNORE`` keeps the earliest (positive *or* negative)
row.
"""

import sqlite3


def pfam_sort_key(pfam):
    """Sort key for a Pfam accession by its numeric part (``PF00028`` -> ``28``).

    Used to canonicalise DDI pairs so the stored column order depends only on the
    Pfam accessions, never on the internal ``domain.id`` insertion order.  Strips
    everything but digits; accessions without any digit fall back to a lexical key
    that sorts deterministically after all numbered ones.
    """
    digits = "".join(c for c in pfam if c.isdigit())
    return (0, int(digits)) if digits else (1, pfam)


def ensure_domains(conn, pfam_ids):
    """Bulk ``INSERT OR IGNORE`` domain rows for ``pfam_ids`` (name left NULL).

    Returns the number of distinct Pfam IDs supplied.
    """
    unique = {p for p in pfam_ids if p}
    conn.executemany(
        "INSERT OR IGNORE INTO domain(pfam_id) VALUES (?)",
        [(p,) for p in unique],
    )
    return len(unique)


def _pfam_to_id(conn):
    return {pfam: did for did, pfam in conn.execute("SELECT id, pfam_id FROM domain")}


def insert_ddis(conn, pairs, negative, source):
    """Insert DDIs for ``(pfam_a, pfam_b)`` pairs.

    Domains must already exist (call :func:`ensure_domains` first); pairs whose
    Pfam is missing from the ``domain`` table are skipped.  Each pair is stored
    as ``(min(id), max(id))`` so swapped duplicates collapse onto one row.

    Returns the number of rows offered to ``INSERT OR IGNORE`` (before dedup by
    the DB).
    """
    pfam_to_id = _pfam_to_id(conn)
    neg = int(bool(negative))
    rows = []
    seen = set()
    for a, b in pairs:
        ia = pfam_to_id.get(a)
        ib = pfam_to_id.get(b)
        if ia is None or ib is None:
            continue
        key = (ia, ib) if pfam_sort_key(a) <= pfam_sort_key(b) else (ib, ia)
        if key in seen:
            continue
        seen.add(key)
        rows.append((key[0], key[1], neg, source))

    conn.executemany(
        "INSERT OR IGNORE INTO domain_domain_interaction"
        "(domain_id_a, domain_id_b, negative, source) VALUES (?, ?, ?, ?)",
        rows,
    )
    return len(rows)


def count_source(conn, source):
    """Number of DDI rows currently tagged with ``source``."""
    return conn.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction WHERE source = ?",
        (source,),
    ).fetchone()[0]
