#!/usr/bin/env python3
"""BUILD_EXTERNAL_TEST must give the same answer either side of enrichment.

`workflows/domainsplit.nf` schedules BUILD_EXTERNAL_TEST *before*
ENRICH_DDI_DATABASE, because it clones the master and after enrichment ~99% of
that file is per-residue ProtT5/ESM blobs this step never reads (590 MB of
594 MB at `-profile test` scale, tens of GB in production).

That reorder is only legitimate if it cannot change the output, and the argument
is narrow enough to be worth testing rather than reasoning about: the two steps
write disjoint tables, and their one shared table is `domain_protein_map`, which
`insert_domain_protein_mapping.py` upserts. It leaves `instance_id` NULL on rows
it adds and updates only `domain_sequence` / `esm*_per_domain` on rows it hits,
while `build_external_test.py` reads `domain_protein_map` WHERE
`instance_id IS NOT NULL`. So the enrichment step is invisible to it.

Asserted here: running build_external_test.py against the pre-enrichment DB and
against the same DB after the two enrichment steps that touch `protein` and
`domain_protein_map` yields byte-identical `ddi_split_membership` rows and an
identical drop report. The enrichment steps are the real scripts, not a
simulation of them -- a change to their upsert that started writing
`instance_id`, or adding instance rows, would break this test, which is exactly
the regression the reorder is exposed to.

Run directly (`python3 tests/python/test_external_test_before_enrich.py`) or via
pytest.
"""

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
INSERT_PROTEINS = os.path.join(BIN, "insert_proteins_with_embeddings.py")
INSERT_MAPPING = os.path.join(BIN, "insert_domain_protein_mapping.py")

from ddi_db_utils import ensure_domains, insert_ddis  # noqa: E402
from test_enrich_after_ingest import INSTANCES, build_inputs  # noqa: E402

ENV = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))

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
        db, pd_map, fasta, prott5, esm_protein, esm_domain = build_inputs(tmp)
        add_external_ddis(db)

        before = run_external_test(db, tmp, "before")

        # The two enrichment steps that touch `protein` and `domain_protein_map`;
        # the GO / PPI steps write tables build_external_test.py never reads.
        # Both copy --db-in to ./domainsplit.sqlite3, so chain them through a
        # staged name the way ENRICH_DDI_DATABASE chains the processes.
        staged = os.path.join(tmp, "staged.sqlite3")
        shutil.copyfile(db, staged)
        for cmd in (
            [sys.executable, INSERT_PROTEINS, "--db-in", staged, "--uniprot-db", fasta,
             "--protein-domain-map", pd_map, "--prott5-embeddings", prott5,
             "--esm-protein-embeddings", esm_protein,
             "--versions", os.path.join(tmp, "v1.yml"), "--process-name", "TEST"],
            [sys.executable, INSERT_MAPPING, "--db-in", staged,
             "--protein-domain-map", pd_map, "--esm-domain-embeddings", esm_domain,
             "--versions", os.path.join(tmp, "v2.yml"), "--process-name", "TEST"],
        ):
            done = subprocess.run(cmd, env=ENV, cwd=tmp, capture_output=True, text=True)
            assert done.returncode == 0, done.stderr
            os.replace(os.path.join(tmp, "domainsplit.sqlite3"), staged)
        enriched = staged

        # Enrichment really did happen -- otherwise this test would pass vacuously.
        conn = sqlite3.connect(enriched)
        n_blobs = conn.execute(
            "SELECT COUNT(*) FROM protein WHERE esm3_per_residue IS NOT NULL"
        ).fetchone()[0]
        n_instances = conn.execute(
            "SELECT COUNT(*) FROM domain_protein_map WHERE instance_id IS NOT NULL"
        ).fetchone()[0]
        conn.close()
        assert n_blobs > 0, "per-residue embeddings were not written"
        assert n_instances == len(INSTANCES), n_instances

        after = run_external_test(enriched, tmp, "after")
        assert after == before, f"\nbefore={before}\nafter ={after}"
        assert before[0], "no membership rows were produced at all"


if __name__ == "__main__":
    test_membership_is_unchanged_by_enrichment()
    print("OK: external test membership is identical either side of enrichment")
