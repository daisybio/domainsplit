#!/usr/bin/env python3
"""Checks for the structural steps: enrich_structural_af.py, score_ddi.py's
sampling and pair selection, and merge_scores.py.

Asserted here:

  * ENRICH_STRUCTURAL slices an instance pair out of its complex with the
    chains in the *row's* orientation, even when af_metadata lists the pair the
    other way round, and sets has_structure = 1 on every row of that triple --
    including rows of a method it was not asked to slice, because the slice is
    in structures.h5 either way. A pair without metadata keeps the default 0.
  * score_ddi.py scores each (ddi, instance pair) of a method once, however
    many splits hold it (the 1:1 and 1:10 test sets share their positives),
    and its random trials depend only on (seed, pair), not on shard or order.
  * merge_scores.py writes that one z-score to every split row of the triple,
    for its own method only.

Run directly (`python3 tests/python/test_structure_scoring.py`) or via pytest.
"""

import gzip
import io
import os
import sqlite3
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from domainsplit_schema import add_instance, make_db  # noqa: E402

ENV = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))
N_RES = 20


def pdb_text(chains):
    """Minimal two-chain PDB: ALA residues 1..N_RES along x, one chain per y offset."""
    lines, serial = [], 1
    for chain_id, y in chains:
        for resseq in range(1, N_RES + 1):
            x = 3.8 * resseq
            for name, dx, dz, element in (("N", -1.2, 0.0, "N"), ("CA", 0.0, 0.0, "C"),
                                          ("C", 1.2, 0.0, "C"), ("O", 1.2, 1.2, "O"),
                                          ("CB", 0.0, -1.5, "C")):
                lines.append(
                    f"ATOM  {serial:5d} {name:<4s} ALA {chain_id}{resseq:4d}    "
                    f"{x + dx:8.3f}{y:8.3f}{dz:8.3f}{1.0:6.2f}{90.0:6.2f}          {element:>2s}"
                )
                serial += 1
        lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"


def build(tmp):
    import zstandard

    db = os.path.join(tmp, "in.sqlite3")
    make_db(db)
    conn = sqlite3.connect(db)
    ia = add_instance(conn, "PF00001", "P1", 5, 15)
    ib = add_instance(conn, "PF00002", "P2", 3, 9)
    ic = add_instance(conn, "PF00003", "P3", 1, 10)
    conn.execute("INSERT INTO domain_domain_interaction(domain_id_a, domain_id_b, negative, source) "
                 "VALUES (1, 2, 0, '3did'), (1, 3, 1, 'sampled_negative')")
    rows = [(1, "random", s, ia, ib) for s in ("train", "test_uniform_1_1", "test_uniform_1_10")]
    rows += [(1, "random_ilp", "train", ia, ib), (2, "random", "train", ia, ic)]
    conn.executemany("INSERT INTO ddi_split_membership(ddi_id, method, split, instance_id_a, "
                     "instance_id_b) VALUES (?, ?, ?, ?, ?)", rows)
    conn.commit()
    conn.close()

    # The file has P2 as chain A and P1 as chain B -- the reverse of the row.
    # P1's chain sits at y=0 and P2's at y=4, so the slice's chain A must be at y=0.
    complex_path = os.path.join(tmp, "p2_p1.pdb.zst")
    with open(complex_path, "wb") as fh:
        fh.write(zstandard.ZstdCompressor().compress(pdb_text([("A", 4.0), ("B", 0.0)]).encode()))
    meta = os.path.join(tmp, "af_metadata.csv")
    with open(meta, "w") as fh:
        fh.write("uniprot_id_a,uniprot_id_b,path,source\n")
        # Relative to the CSV, as a fixture has to write it.
        fh.write(f"P2,P1,{os.path.basename(complex_path)},downloaded\n")
    return db, meta, ia, ib, ic


def run_enrich(tmp):
    db, meta, ia, ib, ic = build(tmp)
    out_db = os.path.join(tmp, "out.sqlite3")
    h5_path = os.path.join(tmp, "structures.h5")
    subprocess.run(
        [sys.executable, os.path.join(BIN, "enrich_structural_af.py"),
         "--db_in", db, "--db_out", out_db, "--structures_h5", h5_path,
         "--af_metadata", meta, "--methods", "random", "--cpus", "1",
         "--versions", os.path.join(tmp, "versions.yml"), "--process_name", "TEST"],
        check=True, env=ENV, capture_output=True, text=True,
    )
    return out_db, h5_path, ia, ib, ic


def test_enrich_slices_in_row_orientation_and_flags_every_row():
    import h5py
    from Bio.PDB.PDBParser import PDBParser

    import utils_struct

    with tempfile.TemporaryDirectory() as tmp:
        out_db, h5_path, ia, ib, ic = run_enrich(tmp)
        conn = sqlite3.connect(out_db)
        flags = dict(((m, s, d), f) for d, m, s, f in conn.execute(
            "SELECT ddi_id, method, split, has_structure FROM ddi_split_membership"))
        conn.close()
        assert flags[("random", "train", 1)] == 1
        assert flags[("random", "test_uniform_1_10", 1)] == 1
        assert flags[("random_ilp", "train", 1)] == 1, "has_structure is per instance pair, not per method"
        assert flags[("random", "train", 2)] == 0, "no metadata -> no structure"

        with h5py.File(h5_path, "r") as h5:
            raw = utils_struct.get_structure_bytes(h5, ia, ib)
            assert raw is not None
            assert utils_struct.get_structure_bytes(h5, ia, ic) is None
        s = PDBParser(QUIET=True).get_structure("x", io.StringIO(gzip.decompress(raw).decode()))
        a, b = s[0]["A"], s[0]["B"]
        assert [r.id[1] for r in a] == list(range(5, 16))
        assert [r.id[1] for r in b] == list(range(3, 10))
        assert a[5]["CA"].coord[1] == 0.0, "chain A must be P1 (y=0) after the swap"
        assert b[3]["CA"].coord[1] == 4.0


def test_scoring_pairs_are_distinct_and_draws_are_seeded():
    import score_ddi
    import utils_struct

    with tempfile.TemporaryDirectory() as tmp:
        out_db, _h5, ia, ib, _ic = run_enrich(tmp)
        rows = utils_struct.get_structures_for_method(out_db, "random")
        assert rows == [(1, ia, ib)], "three split rows of one pair must be scored once"
        assert utils_struct.get_positive_train_structures(out_db, "random") == [(1, ia, ib)]

    draw = lambda seed, a, b: score_ddi.pair_rng(seed, a, b).integers(0, 1 << 30, 5).tolist()  # noqa: E731
    assert draw(42, "x", "y") == draw(42, "x", "y")
    assert draw(42, "x", "y") != draw(43, "x", "y")
    assert draw(42, "x", "y") != draw(42, "y", "x")


def test_merge_scores_copies_z_to_every_split_row_of_its_method():
    with tempfile.TemporaryDirectory() as tmp:
        out_db, _h5, ia, ib, _ic = run_enrich(tmp)
        scores = os.path.join(tmp, "random.scores.tsv")
        confirmed = os.path.join(tmp, "random.confirmed.tsv")
        with open(scores, "w") as fh:
            fh.write("ddi_id\tmethod\tinstance_id_a\tinstance_id_b\tz_score\n")
            fh.write(f"1\trandom\t{ia}\t{ib}\t3.5\n")
        with open(confirmed, "w") as fh:
            fh.write("ddi_id\tmethod\tmajority_confirmed\tmean_confirmed\n1\trandom\t1\t1\n")
        merged = os.path.join(tmp, "merged.sqlite3")
        subprocess.run(
            [sys.executable, os.path.join(BIN, "merge_scores.py"), "--db_in", out_db,
             "--db_out", merged, "--scores", scores, "--confirmed", confirmed,
             "--versions", os.path.join(tmp, "v.yml"), "--process_name", "TEST"],
            check=True, env=ENV, capture_output=True, text=True,
        )
        conn = sqlite3.connect(merged)
        z = dict(((m, s, d), v) for d, m, s, v in conn.execute(
            "SELECT ddi_id, method, split, z_score FROM ddi_split_membership"))
        votes = conn.execute("SELECT * FROM ddi_interaction_confirmed").fetchall()
        conn.close()
        assert {z[("random", s, 1)] for s in ("train", "test_uniform_1_1", "test_uniform_1_10")} == {3.5}
        assert z[("random_ilp", "train", 1)] is None, "another method's row keeps its own (absent) score"
        assert votes == [(1, "random", 1, 1)]


def test_rsa_normalises_with_the_requested_chains_own_residue_type():
    from Bio.PDB import PDBParser
    from utils_struct import calculate_rsa_residue_level

    # Chain B renumbered 1..N_RES as in AF3 complexes, but TRP: before the fix
    # its residues overwrote chain A's and A was divided by TRP's max SASA.
    text = pdb_text([("A", 0.0), ("B", 30.0)])
    text = "\n".join(l.replace(" ALA B", " TRP B") for l in text.splitlines()) + "\n"
    st = PDBParser(QUIET=True).get_structure("x", io.StringIO(text))
    rsa = calculate_rsa_residue_level(st, "A")
    assert len(rsa) == N_RES
    sasa = {r.get_id(): r.sasa for r in st[0]["A"]}
    for rid, v in rsa.items():
        assert abs(v - sasa[rid] / 129.0) < 1e-9, "chain A is ALA: normalise by ALA's 129"


if __name__ == "__main__":
    test_rsa_normalises_with_the_requested_chains_own_residue_type()
    test_enrich_slices_in_row_orientation_and_flags_every_row()
    test_scoring_pairs_are_distinct_and_draws_are_seeded()
    test_merge_scores_copies_z_to_every_split_row_of_its_method()
    print("ok")
