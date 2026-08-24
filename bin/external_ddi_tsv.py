#!/usr/bin/env python3
"""The normalized TSV that every external DDI source is parsed into.

External sources (``single_domain_ppi``, ``PPIDM``/``PPIDM_<Class>``,
``negatome``) are not written to the database when they are parsed: their rows
are inserted last, after ppi-splitting has returned, so protected sources keep
first claim on a pair (see ``bin/ddi_db_utils.py``).  Their Pfam families are
needed early though -- they feed the union domain fetch -- so parsing and
inserting are split, with this file as the interchange format.

Columns: ``pfam_a  pfam_b  negative  source``.  Pairs are written in the
canonical order (:func:`ddi_db_utils.canonical_pair`); ``negative`` is 0 or 1.
One physical row per (pair, source): a source that contributes a pair under two
labels -- PPIDM's umbrella ``PPIDM`` plus its per-class ``PPIDM_Gold`` -- writes
two rows, which the inserter merges into one DB row.
"""

import csv

from ddi_db_utils import canonical_pair

COLUMNS = ["pfam_a", "pfam_b", "negative", "source"]


def write_external_tsv(path, rows):
    """Write ``(pfam_a, pfam_b, negative, source)`` rows; returns the count written.

    Rows whose Pfam accessions are missing are skipped.  Exact duplicates are
    collapsed -- the same pair under the same source and label is one row.
    """
    seen, out = set(), []
    for pfam_a, pfam_b, negative, source in rows:
        if not pfam_a or not pfam_b or not source:
            continue
        a, b = canonical_pair(pfam_a, pfam_b)
        key = (a, b, int(bool(negative)), source)
        if key in seen:
            continue
        seen.add(key)
        out.append(key)

    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(COLUMNS)
        writer.writerows(out)
    return len(out)


def read_external_tsv(path):
    """Yield ``(pfam_a, pfam_b, negative, source)`` tuples from ``path``."""
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = set(COLUMNS) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        for row in reader:
            yield (
                row["pfam_a"].strip(),
                row["pfam_b"].strip(),
                int(row["negative"]),
                row["source"].strip(),
            )
