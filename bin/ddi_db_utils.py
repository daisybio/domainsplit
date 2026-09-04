#!/usr/bin/env python3
"""Shared helpers for inserting DDIs into the domainsplit SQLite.

One row per domain pair
-----------------------
``domain_domain_interaction`` is keyed ``UNIQUE(domain_id_a, domain_id_b)``: a
pair of Pfam families has exactly one row, whose ``source`` is a comma-joined
list of every source that contributed it (``"single_domain_ppi,PPIDM,PPIDM_Gold"``).
Pairs are order-normalised by Pfam accession number (see :func:`pfam_sort_key`)
so the stored ``(domain_id_a, domain_id_b)`` order depends only on the
accessions, never on ``domain.id`` insertion order.

Two classes of source
---------------------
* **Protected** (:data:`PROTECTED_SOURCES`) -- ``3did`` and ``sampled_negative``.
  These are the population that ppi-splitting partitions, so they own their
  pairs outright: an incoming external row for a pair already held by a
  protected source is dropped silently.  This is what makes the external test
  set strictly unseen -- a pair dropped here never appeared in *any* split of
  *any* set.
* **External** -- ``single_domain_ppi``, ``PPIDM``/``PPIDM_<Class>``, ``negatome``.
  These merge with each other, and are inserted last.

Insert rules (:func:`merge_external_ddis`)
------------------------------------------
=============================  ==========================================
incoming vs. existing row      action
=============================  ==========================================
no existing row                insert
existing is protected          drop incoming silently
existing external, same label  append source to the comma list
existing external, other label **drop the pair entirely** and log it
=============================  ==========================================

The last rule also applies within one batch: if two external sources disagree
on ``negative`` for the same pair, neither is kept and any existing external row
for that pair is deleted.  Resolution is therefore independent of the order the
sources are supplied in, and a third source cannot resurrect a poisoned pair.
Dropped contributions are reported one row per contribution so the log is
greppable by source.
"""

import sqlite3

#: Sources that own their pairs outright; external sources never merge into them.
PROTECTED_SOURCES = frozenset({"3did", "sampled_negative"})

#: Column order of the conflict report written by INSERT_EXTERNAL_SOURCES.
CONFLICT_COLUMNS = ["pfam_a", "pfam_b", "source", "negative", "action"]

ACTION_DROPPED = "dropped_label_conflict"
ACTION_DELETED = "deleted_existing_label_conflict"


def pfam_sort_key(pfam):
    """Sort key for a Pfam accession by its numeric part (``PF00028`` -> ``28``).

    Used to canonicalise DDI pairs so the stored column order depends only on the
    Pfam accessions, never on the internal ``domain.id`` insertion order.  Strips
    everything but digits; accessions without any digit fall back to a lexical key
    that sorts deterministically after all numbered ones.
    """
    digits = "".join(c for c in pfam if c.isdigit())
    return (0, int(digits)) if digits else (1, pfam)


def canonical_pair(pfam_a, pfam_b):
    """``(pfam_a, pfam_b)`` ordered by :func:`pfam_sort_key`."""
    return tuple(sorted((pfam_a, pfam_b), key=pfam_sort_key))


def split_sources(source):
    """``"PPIDM,PPIDM_Gold"`` -> ``["PPIDM", "PPIDM_Gold"]`` (empty-safe)."""
    if not source:
        return []
    return [s for s in (part.strip() for part in source.split(",")) if s]


def dedup_sources(sources):
    """``sources`` with duplicates and empties removed, first-seen order preserved."""
    seen, ordered = set(), []
    for s in sources:
        if s and s not in seen:
            seen.add(s)
            ordered.append(s)
    return ordered


def join_sources(sources):
    """Comma-join ``sources``, dropping duplicates and preserving first-seen order."""
    return ",".join(dedup_sources(sources))


def is_protected(source):
    """True if any token of ``source`` is a protected source."""
    return any(s in PROTECTED_SOURCES for s in split_sources(source))


def ensure_domains(conn, pfam_ids):
    """Bulk ``INSERT OR IGNORE`` domain rows for ``pfam_ids`` (name left NULL).

    Returns the number of distinct Pfam IDs supplied.  Domains created here for a
    family that never receives an instance are removed later by the prune step,
    together with the DDIs referencing them.
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
    """Insert DDIs for a protected source (``3did``, ``sampled_negative``).

    Domains must already exist (call :func:`ensure_domains` first); pairs whose
    Pfam is missing from the ``domain`` table are skipped.  Each pair is stored
    once, canonicalised; a pair already present -- from this source or an earlier
    protected one -- is left untouched (``INSERT OR IGNORE``).

    Returns the number of rows offered to the DB (before its own dedup).
    """
    pfam_to_id = _pfam_to_id(conn)
    neg = int(bool(negative))
    rows, seen = [], set()
    for a, b in pairs:
        ia, ib = pfam_to_id.get(a), pfam_to_id.get(b)
        if ia is None or ib is None:
            continue
        first, second = canonical_pair(a, b)
        key = (pfam_to_id[first], pfam_to_id[second])
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


def _existing_rows(conn):
    """``{(pfam_a, pfam_b): (ddi_id, negative, source)}`` for every stored DDI."""
    out = {}
    for ddi_id, pfam_a, pfam_b, negative, source in conn.execute(
        "SELECT ddi.id, da.pfam_id, db.pfam_id, ddi.negative, ddi.source "
        "FROM domain_domain_interaction AS ddi "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b"
    ):
        out[canonical_pair(pfam_a, pfam_b)] = (ddi_id, negative, source)
    return out


def merge_external_ddis(conn, rows):
    """Apply the external-source insert rules for ``(pfam_a, pfam_b, negative, source)``.

    ``rows`` is consumed in full before anything is written, so the outcome does
    not depend on the order sources are supplied in.  Domains are bulk-created
    for every referenced Pfam.  Returns ``(stats, conflicts)`` where ``stats`` is
    a dict of counters and ``conflicts`` is a list of :data:`CONFLICT_COLUMNS`
    tuples describing every dropped contribution.
    """
    # pair -> {negative: [source, ...]}, first-seen order preserved throughout
    incoming, pfams = {}, set()
    for pfam_a, pfam_b, negative, source in rows:
        if not pfam_a or not pfam_b or not source:
            continue
        pair = canonical_pair(pfam_a, pfam_b)
        pfams.update(pair)
        incoming.setdefault(pair, {}).setdefault(int(bool(negative)), []).append(source)

    ensure_domains(conn, pfams)
    pfam_to_id = _pfam_to_id(conn)
    existing = _existing_rows(conn)

    stats = {
        "pairs_offered": len(incoming),
        "inserted": 0,
        "merged": 0,
        "dropped_protected": 0,
        "conflict_pairs": 0,
        "sources_dropped": 0,
    }
    conflicts = []
    to_insert, to_update, to_delete = [], [], []

    for pair, by_label in incoming.items():
        pfam_a, pfam_b = pair
        ia, ib = pfam_to_id.get(pfam_a), pfam_to_id.get(pfam_b)
        if ia is None or ib is None:
            continue
        prior = existing.get(pair)

        # A protected source owns the pair; every incoming row is dropped without
        # a report line -- this is the ordinary, expected case, not an anomaly.
        if prior is not None and is_protected(prior[2]):
            stats["dropped_protected"] += 1
            continue

        labels = sorted(by_label)
        conflict = len(labels) > 1 or (prior is not None and prior[1] not in labels)

        if conflict:
            stats["conflict_pairs"] += 1
            for label in labels:
                for source in dedup_sources(by_label[label]):
                    conflicts.append((pfam_a, pfam_b, source, label, ACTION_DROPPED))
                    stats["sources_dropped"] += 1
            if prior is not None:
                to_delete.append((prior[0],))
                for source in split_sources(prior[2]):
                    conflicts.append((pfam_a, pfam_b, source, prior[1], ACTION_DELETED))
            continue

        label = labels[0]
        sources = by_label[label]
        if prior is None:
            to_insert.append((ia, ib, label, join_sources(sources)))
            stats["inserted"] += 1
        else:
            merged = join_sources(split_sources(prior[2]) + sources)
            if merged != prior[2]:
                to_update.append((merged, prior[0]))
            stats["merged"] += 1

    if to_delete:
        conn.executemany(
            "DELETE FROM domain_domain_interaction WHERE id = ?", to_delete
        )
    if to_insert:
        conn.executemany(
            "INSERT OR IGNORE INTO domain_domain_interaction"
            "(domain_id_a, domain_id_b, negative, source) VALUES (?, ?, ?, ?)",
            to_insert,
        )
    if to_update:
        conn.executemany(
            "UPDATE domain_domain_interaction SET source = ? WHERE id = ?", to_update
        )
    return stats, conflicts


def count_source(conn, source):
    """Number of DDI rows whose comma-joined ``source`` list contains ``source``."""
    return conn.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction "
        "WHERE ',' || source || ',' LIKE '%,' || ? || ',%'",
        (source,),
    ).fetchone()[0]


def count_ddis(conn):
    """Total number of DDI rows."""
    return conn.execute(
        "SELECT COUNT(*) FROM domain_domain_interaction"
    ).fetchone()[0]
