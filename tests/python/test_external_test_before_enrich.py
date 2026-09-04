#!/usr/bin/env python3
"""BUILD_EXTERNAL_TEST must give the same answer either side of enrichment.

`workflows/domainsplit.nf` schedules BUILD_EXTERNAL_TEST *before*
ENRICH_DDI_DATABASE, because it clones the master and enrichment only ever grows
it -- so cloning first is the smaller copy. (The margin used to be enormous:
enrichment wrote per-residue embedding blobs that were ~99% of the file. Those are
gone, published as HDF5 instead, and the ordering is now merely the cheaper one.)

That reorder is only legitimate if it cannot change the output, and the argument
is narrow enough to be worth testing rather than reasoning about:
`build_external_test.py` reads `domain_domain_interaction`, `domain` and
`domain_protein_map` WHERE `instance_id IS NOT NULL`, and enrichment now writes
none of those -- INSERT_DOMAIN_PROTEIN_MAPPING went with the embedding columns, so
`domain_protein_map` is written once, by ingest, and never touched again.

Asserted here: running build_external_test.py against the pre-enrichment DB and
against the same DB after the enrichment step that touches `protein` yields
byte-identical `ddi_split_membership` rows and an identical drop report. The
enrichment step is the real script, not a simulation of it -- a change that
started adding `domain_protein_map` rows would break this test, which is exactly
the regression the reorder is exposed to.

Run directly (`python3 tests/python/test_external_test_before_enrich.py`) or via
pytest.
"""

import gzip
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

BUILD_EXTERNAL = os.path.join(BIN, "build_external_test.py")
INSERT_PROTEINS = os.path.join(BIN, "insert_protein_sequences.py")

from ddi_db_utils import ensure_domains, insert_ddis  # noqa: E402
from domainsplit_schema import add_instance, make_db  # noqa: E402

ENV = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))

# Deliberately *not* imported from test_enrich_after_ingest: that file's fixture
# exists to exercise the protein upsert, and coupling the two made a change to one
# test silently reshape the other.
INSTANCES = [
    ("PF00001", "P11111", 10, 20),
    ("PF00001", "P11111", 50, 60),
    ("PF00002", "Q22222", 5, 15),
]
SEQUENCES = {"P11111": "M" * 80, "Q22222": "K" * 40}


def build_inputs(tmp):
    """`(db, protein_domain_map, uniprot_fasta)` for the instances above."""
    db = os.path.join(tmp, "domainsplit.sqlite3")
    make_db(db)
    conn = sqlite3.connect(db)
    conn.execute("PRAGMA foreign_keys=ON")
    for pfam, uniprot, start, end in INSTANCES:
        add_instance(conn, pfam, uniprot, start, end, clan="CL0192", sequence="A" * (end - start + 1))
    conn.commit()
    conn.close()

    pd_map = os.path.join(tmp, "protein_domain_mapping.csv.gz")
    with gzip.open(pd_map, "wt", newline="") as fh:
        fh.write("pfam_id,uniprot_id,start_pos,end_pos,sequence\n")
        for pfam, uniprot, start, end in INSTANCES:
            fh.write(f"{pfam},{uniprot},{start},{end},{'A' * (end - start + 1)}\n")

    fasta = os.path.join(tmp, "uniprot_sequences.fasta.gz")
    with gzip.open(fasta, "wt") as fh:
        for uniprot, seq in SEQUENCES.items():
            fh.write(f">sp|{uniprot}|{uniprot}_HUMAN test protein\n{seq}\n")

    return db, pd_map, fasta

# Non-protected sources, so build_external_test.py treats them as the held-out
# set. Both families carry instances via INSTANCES above.
EXTERNAL_PAIRS = [("PF00001", "PF00002"), ("PF00001", "PF00001")]


def add_external_ddis(db):
    conn = sqlite3.connect(db)
    conn.execute("PRAGMA foreign_keys=ON")
    ensure_domains(conn, [p for pair in EXTERNAL_PAIRS for p in pair])
    insert_ddis(conn, EXTERNAL_PAIRS, negative=False, source="PPIDM_Gold")
    conn.commit()
    conn.close()


def run_external_test(db, tmp, tag):
    out_db = os.path.join(tmp, f"{tag}.sqlite3")
    shutil.copyfile(db, out_db)
    dropped = os.path.join(tmp, f"{tag}_dropped.tsv")
    done = subprocess.run(
        [sys.executable, BUILD_EXTERNAL, "--db", out_db, "--method", "external_test",
         "--split", "test", "--target", "2", "--pool-factor", "3", "--seed", "42",
         "--dropped-out", dropped, "--versions", os.path.join(tmp, f"{tag}_versions.yml"),
         "--process-name", "TEST"],
        env=ENV, capture_output=True, text=True,
    )
    assert done.returncode == 0, done.stderr

    conn = sqlite3.connect(out_db)
    # Resolved to accessions and instance ids, not raw row ids: the point is that
    # the *same* instances are chosen, and enrichment may renumber nothing but
    # comparing surrogate keys would hide a real difference behind a match.
    membership = sorted(conn.execute(
        "SELECT da.pfam_id, db.pfam_id, m.method, m.split, m.instance_id_a, m.instance_id_b "
        "FROM ddi_split_membership AS m "
        "JOIN domain_domain_interaction AS ddi ON ddi.id = m.ddi_id "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b"
    ).fetchall())
    conn.close()
    return membership, open(dropped).read()


def test_membership_is_unchanged_by_enrichment():
    with tempfile.TemporaryDirectory() as tmp:
        db, pd_map, fasta = build_inputs(tmp)
        add_external_ddis(db)

        before = run_external_test(db, tmp, "before")

        # The one enrichment step that touches a table build_external_test.py
        # reads a *neighbour* of; the GO / PPI steps write tables it never reads.
        # It copies --db-in to ./domainsplit.sqlite3, so stage it under another
        # name the way ENRICH_DDI_DATABASE chains the processes.
        staged = os.path.join(tmp, "staged.sqlite3")
        shutil.copyfile(db, staged)
        done = subprocess.run(
            [sys.executable, INSERT_PROTEINS, "--db-in", staged, "--uniprot-db", fasta,
             "--protein-domain-map", pd_map,
             "--versions", os.path.join(tmp, "v1.yml"), "--process-name", "TEST"],
            env=ENV, cwd=tmp, capture_output=True, text=True,
        )
        assert done.returncode == 0, done.stderr
        os.replace(os.path.join(tmp, "domainsplit.sqlite3"), staged)
        enriched = staged

        # Enrichment really did happen -- otherwise this test would pass
        # vacuously. There are no blobs left to look for, so the proxy is the
        # column the step owns: `protein.sequence`, NULL until it runs.
        conn = sqlite3.connect(enriched)
        n_sequences = conn.execute(
            "SELECT COUNT(*) FROM protein WHERE sequence IS NOT NULL"
        ).fetchone()[0]
        n_instances = conn.execute(
            "SELECT COUNT(*) FROM domain_protein_map WHERE instance_id IS NOT NULL"
        ).fetchone()[0]
        conn.close()
        assert n_sequences == len(SEQUENCES), n_sequences
        assert n_instances == len(INSTANCES), n_instances

        after = run_external_test(enriched, tmp, "after")
        assert after == before, f"\nbefore={before}\nafter ={after}"
        assert before[0], "no membership rows were produced at all"


if __name__ == "__main__":
    test_membership_is_unchanged_by_enrichment()
    print("OK: external test membership is identical either side of enrichment")
