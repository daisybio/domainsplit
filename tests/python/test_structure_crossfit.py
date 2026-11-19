#!/usr/bin/env python3
"""Checks for score_ddi.py's potential handling:

  * a sparse potential (no cell reaches E_ab >= 5) gives constant random
    trials, so the pair is left unscored (z NULL, no vote) instead of z = 0.0;
  * cross-fit: the potential that scores a train row has its own fold's
    contacts removed, so a train positive never scores against itself; a
    test pair of the same DDI still uses the full potential.
"""

import json
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from test_score_ddi_shards import BIN, chains, potential, run, write_fixture  # noqa: E402

sys.path.insert(0, BIN)


def zero_folds(tmp):
    path = os.path.join(tmp, "folds.json")
    with open(path, "w") as fh:
        json.dump({"c_ab": np.zeros((5, 20, 20)).tolist(), "surface": np.zeros((5, 20)).tolist()}, fh)
    return path


def score_rows(tmp, db, h5, mats, folds):
    s, c = os.path.join(tmp, "s.tsv"), os.path.join(tmp, "c.tsv")
    run("score_ddi.py", "--db_in", db, "--structures_h5", h5, "--method", "m",
        "--c_ab_matrix", mats[0], "--db_freq", mats[1], "--t_db", mats[2], "--fold_counts", folds,
        "--scores_out", s, "--confirmed_out", c, "--seed", "7",
        "--versions", os.path.join(tmp, "v.yml"), "--process_name", "T")
    return open(s).read().splitlines()[1:], open(c).read().splitlines()[1:]


def test_constant_random_trials_leave_the_pair_unscored():
    with tempfile.TemporaryDirectory() as tmp:
        db, h5 = write_fixture(tmp, {1: 0}, [(1, "test", "a1", "b1", chains("TRP"))])
        folds = zero_folds(tmp)
        # T_DB 400, uniform freq: E_ab = 1 < 5 everywhere, so only fallbacks remain.
        assert score_rows(tmp, db, h5, potential(tmp, 400), folds) == ([], [])
        scores, votes = score_rows(tmp, db, h5, potential(tmp, 2000), folds)  # control: E_ab = 5
        assert len(scores) == 1 and len(votes) == 1


def test_train_row_potential_excludes_its_own_folds_contacts():
    import score_ddi
    import utils_struct

    seed = 7
    assert utils_struct.ddi_fold(seed, 1) != utils_struct.ddi_fold(seed, 2)
    with tempfile.TemporaryDirectory() as tmp:
        # DDI 1's only contacts are TRP-TRP, DDI 2's ALA-ALA; both train positives.
        db, h5 = write_fixture(tmp, {1: 0, 2: 0}, [
            (1, "train", "a1", "b1", chains("TRP")), (2, "train", "a2", "b2", chains("ALA"))])
        c, s = os.path.join(tmp, "c.csv"), os.path.join(tmp, "s.csv")
        run("build_scoring_matrix.py", "--db_in", db, "--structures_h5", h5, "--method", "m",
            "--c_ab_counts_out", c, "--surface_counts_out", s, "--seed", str(seed),
            "--versions", os.path.join(tmp, "v.yml"), "--process_name", "T")
        out = {k: os.path.join(tmp, k) for k in ("cab.csv", "freq.csv", "t_db.txt", "folds.json")}
        run("merge_scoring_matrix.py", "--method", "m", "--c_ab_counts", c, "--surface_counts", s,
            "--c_ab_matrix_out", out["cab.csv"], "--db_freq_out", out["freq.csv"],
            "--t_db_out", out["t_db.txt"], "--fold_counts_out", out["folds.json"])

        c_ab, freq = score_ddi.load_c_ab_matrix(out["cab.csv"]), score_ddi.load_db_freq(out["freq.csv"])
        t_db = score_ddi.load_t_db(out["t_db.txt"])
        c_folds, s_folds = score_ddi.load_fold_counts(out["folds.json"])
        trp, ala = (utils_struct.AA_3.index(a) for a in ("TRP", "ALA"))
        assert c_ab[trp, trp] > 0 and c_ab[ala, ala] > 0

        c_1, freq_1, t_1 = score_ddi.holdout_counts(c_ab, freq, t_db, c_folds, s_folds,
                                                    utils_struct.ddi_fold(seed, 1))
        assert c_1[trp, trp] == 0 and freq_1[trp] == 0, "DDI 1 scored against its own contacts"
        assert c_1[ala, ala] == c_ab[ala, ala] and t_1 == c_ab[ala, ala]
        c_2, _, _ = score_ddi.holdout_counts(c_ab, freq, t_db, c_folds, s_folds,
                                             utils_struct.ddi_fold(seed, 2))
        assert c_2[ala, ala] == 0 and c_2[trp, trp] == c_ab[trp, trp]


def test_train_pair_is_cross_fit_and_its_test_twin_uses_the_full_potential():
    import utils_struct

    seed = 7
    k1, k2 = utils_struct.ddi_fold(seed, 1), utils_struct.ddi_fold(seed, 2)
    assert k1 != k2
    with tempfile.TemporaryDirectory() as tmp:
        # DDI 1 has a train pair and a test pair with identical TRP contacts.
        db, h5 = write_fixture(tmp, {1: 0, 2: 0}, [
            (1, "train", "a1", "b1", chains("TRP")), (1, "test", "a2", "b2", chains("TRP")),
            (2, "train", "a3", "b3", chains("TRP"))])
        mats = potential(tmp, 2000)
        # Every TRP-TRP count sits in DDI 1's fold: without it only fallbacks
        # remain, so a pair scored with that potential has constant trials.
        c_folds = np.zeros((5, 20, 20))
        trp = utils_struct.AA_3.index("TRP")
        c_folds[k1, trp, trp] = 1000
        with open(mats[3], "w") as fh:
            json.dump({"c_ab": c_folds.tolist(), "surface": np.zeros((5, 20)).tolist()}, fh)
        scores, _ = score_rows(tmp, db, h5, mats, mats[3])
        scored = {tuple(r.split("\t")[:3]) for r in scores}
        assert ("1", "m", "PF00001_A1_1_6") not in scored, "train pair scored with its own fold's contacts"
        assert ("1", "m", "PF00001_A2_1_6") in scored, "test pair of a train DDI must use the full potential"
        assert ("2", "m", "PF00001_A3_1_6") in scored, "control: another fold's train pair keeps TRP-TRP"


def test_identical_trials_with_rounding_noise_leave_the_pair_unscored():
    import score_ddi

    trials = np.full(1000, 0.1)
    assert np.std(trials) != 0, "fixture no longer exercises float rounding"
    assert score_ddi.compute_z_score(0.5, trials) is None
    assert score_ddi.compute_z_score(0.5, np.array([0.0, 1.0])) == 0.0
