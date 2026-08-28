"""The subset of INIT_DOMAINSPLIT_DB's schema the DDI unit tests need.

Kept in step with ``modules/local/init_domainsplit_db/main.nf`` -- in particular:

* ``UNIQUE(domain_id_a, domain_id_b)`` on ``domain_domain_interaction``, without
  ``source``, which is what makes a domain pair a single row carrying a
  comma-joined source list;
* ``domain_protein_map`` keyed by ``(domain_id, protein_id, start_pos, end_pos)``
  and carrying ``instance_id``, so a protein with two copies of one family is two
  rows -- the embedding H5 key contract is instance-level, and so is
  ppi-splitting. ``start_pos``/``end_pos`` are **INTEGER**: with no declared type
  they take BLOB affinity and TEXT ``'10'`` stops matching INTEGER ``10`` in the
  UNIQUE key, which is exactly the desync this mirror has to reproduce;
* ``ddi_split_membership``, which SUBSET_SPLIT_DB filters on.

The GO / PPI tables are here only so the prune and subset tests can assert that
their cascades fire.
"""

import sqlite3

SCHEMA = """
PRAGMA foreign_keys=ON;

CREATE TABLE domain (id INTEGER PRIMARY KEY, pfam_id, name, UNIQUE(pfam_id));
CREATE TABLE domain_go_terms(
    domain_id REFERENCES domain ON DELETE CASCADE,
    go_accession
);
CREATE TABLE domain_domain_interaction (
    id INTEGER PRIMARY KEY,
    domain_id_a, domain_id_b, negative,
    source VARCHAR(255),
    FOREIGN KEY(domain_id_a) REFERENCES domain ON DELETE CASCADE,
    FOREIGN KEY(domain_id_b) REFERENCES domain ON DELETE CASCADE,
    UNIQUE(domain_id_a, domain_id_b)
);

CREATE TABLE protein (
    id INTEGER PRIMARY KEY,
    uniprot_id,
    sequence,
    UNIQUE(uniprot_id)
);
CREATE TABLE protein_go_terms(
    protein_id REFERENCES protein ON DELETE CASCADE,
    go_accession
);
CREATE TABLE protein_protein_interaction (
    protein_id_a REFERENCES protein ON DELETE CASCADE,
    protein_id_b REFERENCES protein ON DELETE CASCADE,
    score,
    UNIQUE(protein_id_a, protein_id_b)
);

CREATE TABLE domain_protein_map (
    domain_id REFERENCES domain ON DELETE CASCADE,
    protein_id REFERENCES protein ON DELETE CASCADE,
    domain_sequence, start_pos INTEGER, end_pos INTEGER,
    instance_id, clan, taxon_id,
    UNIQUE(domain_id, protein_id, start_pos, end_pos)
);

CREATE TABLE ddi_split_membership (
    ddi_id REFERENCES domain_domain_interaction ON DELETE CASCADE,
    method, split,
    instance_id_a, instance_id_b,
    UNIQUE(ddi_id, method, split, instance_id_a, instance_id_b)
);

CREATE UNIQUE INDEX idx_domain_protein_map_instance_id
ON domain_protein_map (instance_id) WHERE instance_id IS NOT NULL;
"""


def make_db(path):
    """Create an empty domainsplit SQLite at ``path`` and return nothing."""
    conn = sqlite3.connect(path)
    conn.executescript(SCHEMA)
    conn.commit()
    conn.close()


def ddi_rows(conn):
    """``{(pfam_a, pfam_b): (negative, source)}`` for every stored DDI."""
    return {
        (a, b): (negative, source)
        for a, b, negative, source in conn.execute(
            "SELECT da.pfam_id, db.pfam_id, ddi.negative, ddi.source "
            "FROM domain_domain_interaction AS ddi "
            "JOIN domain AS da ON da.id = ddi.domain_id_a "
            "JOIN domain AS db ON db.id = ddi.domain_id_b"
        )
    }


def add_instance(conn, pfam, uniprot, start, end, taxon="9606", clan=None, sequence="AAA"):
    """Insert one domain instance (and its domain/protein rows) and return its instance_id.

    Mirrors what INGEST_INSTANCES writes, including the
    ``{family}_{uniprot}_{start}_{end}`` instance id ppi-splitting uses.

    ``start``/``end`` are bound as **strings** on purpose, so every fixture DB
    exercises the INTEGER affinity coercion instead of hiding it: every caller
    passes Python ints, which is why a fixture DB used to match a TEXT-binding
    writer that the real pipeline did not.
    """
    conn.execute("INSERT OR IGNORE INTO domain(pfam_id) VALUES (?)", (pfam,))
    conn.execute("INSERT OR IGNORE INTO protein(uniprot_id) VALUES (?)", (uniprot,))
    domain_id = conn.execute("SELECT id FROM domain WHERE pfam_id = ?", (pfam,)).fetchone()[0]
    protein_id = conn.execute("SELECT id FROM protein WHERE uniprot_id = ?", (uniprot,)).fetchone()[0]
    instance_id = f"{pfam}_{uniprot}_{start}_{end}"
    conn.execute(
        "INSERT OR IGNORE INTO domain_protein_map"
        "(domain_id, protein_id, domain_sequence, start_pos, end_pos, instance_id, clan, taxon_id) "
        "VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        (domain_id, protein_id, sequence, str(start), str(end), instance_id, clan, taxon),
    )
    return instance_id
