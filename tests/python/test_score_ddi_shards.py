#!/usr/bin/env python3
"""Checks for the structural scoring maths on tiny synthetic structures:

  * build_scoring_matrix.py + merge_scoring_matrix.py give the exact T_DB of a
    known-contact fixture, from this method's train positives only.
  * score_ddi.py gives the same sorted score and confirmation rows with
    --n_shards 1 and 3, on DDIs whose instance pairs vote differently.
  * score_ddi.py with T_DB = 0 writes empty outputs (z_score stays NULL) and
    warns, and merge_scores.py then warns for that method.
"""

import gzip
import os
import sqlite3
import subprocess
import sys
import tempfile

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from domainsplit_schema import make_db  # noqa: E402

ENV = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))


def pdb_text(residues):
    """residues: (chain, resseq, resname, x, y); five backbone+CB atoms each."""
    lines, serial = [], 1
    for chain, resseq, resname, x, y in residues:
        for name, dx, dz, element in (("N", -1.2, 0.0, "N"), ("CA", 0.0, 0.0, "C"),
                                      ("C", 1.2, 0.0, "C"), ("O", 1.2, 1.2, "O"),
                                      ("CB", 0.0, -1.5, "C")):
            lines.append(
                f"ATOM  {serial:5d} {name:<4s} {resname} {chain}{resseq:4d}    "
                f"{x + dx:8.3f}{y:8.3f}{dz:8.3f}{1.0:6.2f}{90.0:6.2f}          {element:>2s}"
            )
            serial += 1
    lines.append("END")
    return "\n".join(lines) + "\n"


def write_fixture(tmp, ddis, pairs):
    """ddis: {ddi_id: negative}; pairs: [(ddi_id, split, ia, ib, pdb)], ia/ib short tags."""
    import h5py

    import utils_struct

    pairs = [(d, s, f"PF00001_{ia.upper()}_1_6", f"PF00002_{ib.upper()}_1_6", pdb) for d, s, ia, ib, pdb in pairs]

    db = os.path.join(tmp, "db.sqlite3")
    make_db(db)
    conn = sqlite3.connect(db)
    conn.executemany("INSERT INTO domain(id, pfam_id) VALUES (?, ?)",
                     [(i, f"PF{i:05d}") for i in range(1, 2 * len(ddis) + 1)])
    conn.executemany("INSERT INTO domain_domain_interaction(id, domain_id_a, domain_id_b, negative, source) "
                     "VALUES (?, ?, ?, ?, 'x')",
                     [(d, 2 * d - 1, 2 * d, neg) for d, neg in ddis.items()])
    conn.executemany("INSERT INTO ddi_split_membership(ddi_id, method, split, instance_id_a, instance_id_b, "
                     "has_structure) VALUES (?, 'm', ?, ?, ?, 1)",
                     [(d, s, ia, ib) for d, s, ia, ib, _ in pairs])
    conn.commit()
    conn.close()
    h5_path = os.path.join(tmp, "structures.h5")
    with h5py.File(h5_path, "w") as h5:
        for _, _, ia, ib, pdb in pairs:
            group_key, dset_key = utils_struct.structure_h5_keys(ia, ib)
            h5.require_group(group_key).create_dataset(
                dset_key, data=np.frombuffer(gzip.compress(pdb.encode()), dtype="uint8"))
    return db, h5_path


def run(script, *args):
    return subprocess.run([sys.executable, os.path.join(BIN, script), *args],
                          check=True, env=ENV, capture_output=True, text=True)


def test_t_db_counts_contacting_residue_pairs_of_train_positives_only():
    # One ALA on chain A; chain B has a TRP 4 Å away (contact) and one 40 Å away.
    contact = pdb_text([("A", 1, "ALA", 0.0, 0.0), ("B", 1, "TRP", 0.0, 4.0), ("B", 2, "TRP", 0.0, 40.0)])
    with tempfile.TemporaryDirectory() as tmp:
        db, h5 = write_fixture(tmp, {1: 0, 2: 0, 3: 1}, [
            (1, "train", "a1", "b1", contact),
            (2, "train", "a2", "b2", contact),
            (2, "test", "a3", "b3", contact),   # not train
            (3, "train", "a4", "b4", contact),  # negative DDI
        ])
        shards = []
        for k in range(2):
            c, s = os.path.join(tmp, f"c{k}.csv"), os.path.join(tmp, f"s{k}.csv")
            run("build_scoring_matrix.py", "--db_in", db, "--structures_h5", h5, "--method", "m",
                "--c_ab_counts_out", c, "--surface_counts_out", s, "--shard_id", str(k), "--n_shards", "2",
                "--versions", os.path.join(tmp, "v.yml"), "--process_name", "T")
            shards.append((c, s))
        t_db = os.path.join(tmp, "t_db.txt")
        run("merge_scoring_matrix.py", "--method", "m", "--c_ab_counts", *[c for c, _ in shards],
            "--surface_counts", *[s for _, s in shards], "--c_ab_matrix_out", os.path.join(tmp, "cab.csv"),
            "--db_freq_out", os.path.join(tmp, "freq.csv"), "--t_db_out", t_db)
        # 2 train positives x 1 contacting (ALA, TRP) pair, counted as (A,B) and (B,A).
        assert open(t_db).read().strip() == "4"


def potential(tmp, t_db):
    from utils_struct import AA_3

    cab, freq, t = (os.path.join(tmp, n) for n in ("cab.csv", "freq.csv", "t_db.txt"))
    with open(cab, "w") as fh:  # only TRP-TRP contacts are favourable
        fh.write("aaA,aaB,count\n" + "".join(
            f"{a},{b},{1000 if a == b == 'TRP' else 0}\n" for a in AA_3 for b in AA_3))
    with open(freq, "w") as fh:
        fh.write("aa,frequency\n" + "".join(f"{a},{1 / len(AA_3)}\n" for a in AA_3))
    with open(t, "w") as fh:
        fh.write(f"{t_db}\n")
    return cab, freq, t


def chains(resname):
    return pdb_text([(c, i, resname, 3.8 * i, y) for c, y in (("A", 0.0), ("B", 4.0)) for i in range(1, 7)])


def score(tmp, db, h5, mats, n_shards):
    scores, confirmed, stderr = [], [], ""
    for k in range(n_shards):
        s, c = os.path.join(tmp, f"{n_shards}_{k}.scores.tsv"), os.path.join(tmp, f"{n_shards}_{k}.confirmed.tsv")
        res = run("score_ddi.py", "--db_in", db, "--structures_h5", h5, "--method", "m",
                  "--c_ab_matrix", mats[0], "--db_freq", mats[1], "--t_db", mats[2],
                  "--scores_out", s, "--confirmed_out", c, "--shard_id", str(k), "--n_shards", str(n_shards),
                  "--seed", "7", "--versions", os.path.join(tmp, "v.yml"), "--process_name", "T")
        scores += open(s).read().splitlines()[1:]
        confirmed += open(c).read().splitlines()[1:]
        stderr += res.stderr
    return sorted(scores), sorted(confirmed), stderr


def test_score_rows_and_votes_do_not_depend_on_shard_count():
    hi, lo = chains("TRP"), chains("ALA")
    with tempfile.TemporaryDirectory() as tmp:
        db, h5 = write_fixture(tmp, {1: 0, 2: 0, 3: 1}, [
            (1, "train", "a1", "b1", hi), (1, "test", "a2", "b2", lo),
            (2, "train", "a3", "b3", hi), (2, "train", "a4", "b4", hi), (2, "val", "a5", "b5", lo),
            (3, "test", "a6", "b6", lo), (3, "test", "a7", "b7", lo),
        ])
        mats = potential(tmp, 2000)
        one, three = score(tmp, db, h5, mats, 1), score(tmp, db, h5, mats, 3)
        assert one[:2] == three[:2]
        assert len(one[0]) == 7
        # DDI 1 splits its votes 1:1: no majority, but the high z lifts the mean.
        assert one[1] == ["1\tm\t0\t1", "2\tm\t1\t1", "3\tm\t0\t0"]


def test_t_db_zero_leaves_scores_null_and_merge_warns():
    with tempfile.TemporaryDirectory() as tmp:
        db, h5 = write_fixture(tmp, {1: 0}, [(1, "train", "a1", "b1", chains("TRP"))])
        scores, confirmed, stderr = score(tmp, db, h5, potential(tmp, 0), 1)
        assert scores == [] and confirmed == []
        assert "method=m has T_DB=0" in stderr
        s, c = os.path.join(tmp, "m.scores.tsv"), os.path.join(tmp, "m.confirmed.tsv")
        os.rename(os.path.join(tmp, "1_0.scores.tsv"), s)
        os.rename(os.path.join(tmp, "1_0.confirmed.tsv"), c)
        merged = os.path.join(tmp, "merged.sqlite3")
        res = run("merge_scores.py", "--db_in", db, "--db_out", merged, "--scores", s, "--confirmed", c,
                  "--versions", os.path.join(tmp, "v.yml"), "--process_name", "T")
        assert "method=m has 1 instance pairs with structure but none scored" in res.stderr
        conn = sqlite3.connect(merged)
        assert conn.execute("SELECT z_score FROM ddi_split_membership").fetchall() == [(None,)]
        conn.close()
