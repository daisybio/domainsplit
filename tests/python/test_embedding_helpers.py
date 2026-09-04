#!/usr/bin/env python3
"""The two things the ProtT5 branch is easy to get silently wrong.

Neither needs a GPU, model weights, or torch: `bin/run_embeddings.py` keeps its
heavy imports inside functions, so the helpers below import on their own. There
is no marker, skip or GPU-fake mechanism anywhere in `tests/` -- rather than add
one, this test only touches the parts that do not need it.

Both failures are silent. Feed ProtT5 a bare, unsubstituted sequence and it
tokenizes to something the checkpoint was never trained on; average the `</s>`
position into the pooled vector and every embedding is contaminated by a
separator token. Neither raises, and both would only show up as a model that
learns slightly less well than it should.

Run directly (`python3 tests/python/test_embedding_helpers.py`) or via pytest.
"""

import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "bin"))

from run_embeddings import prepare_prott5_sequence, real_token_lengths  # noqa: E402


def test_rare_residues_become_x():
    # ProtT5's own recipe: U, Z, O and B are mapped to X. Nothing else changes.
    assert prepare_prott5_sequence("ACDUZOBEFG") == "A C D X X X X E F G"
    assert prepare_prott5_sequence("MKV") == "M K V"
    # Lower case reaches this from a FASTA now and then; the vocabulary is upper.
    assert prepare_prott5_sequence("mkvUx") == "M K V X X"


def test_residues_are_whitespace_separated():
    seq = "MKVLIA"
    prepared = prepare_prott5_sequence(seq)
    assert prepared.split(" ") == list(seq)
    # One token per residue and no leading/trailing space, or the tokenizer emits
    # a different number of positions than the mask arithmetic below assumes.
    assert len(prepared.split(" ")) == len(seq)
    assert prepared == prepared.strip()


def test_eos_is_excluded_from_the_pooled_mean():
    # T5 appends exactly one `</s>`, so attention_mask.sum(1) is residues + 1 and
    # the mean must run over `sum - 1` positions.
    assert real_token_lengths([7]) == [6]
    assert real_token_lengths([31, 11, 4]) == [30, 10, 3]
    # A tensor-derived sum arrives as a float or a numpy scalar depending on the
    # torch version; the helper takes int() itself.
    assert real_token_lengths([7.0]) == [6]


def test_pooled_length_never_reaches_zero():
    # An all-padding row would otherwise divide by zero in masked_mean.
    assert real_token_lengths([1]) == [1]
    assert real_token_lengths([0]) == [1]


def test_esm_and_prott5_pool_different_spans_on_purpose():
    """The asymmetry is deliberate and worth pinning down.

    The ESM branch pools over the full tokenized length (BOS + seq + EOS); ProtT5
    pools over residues only. If someone ever makes them agree, one of the two
    stops matching its model's own embedding recipe -- so assert that a 6-residue
    sequence gives 6 ProtT5 positions while ESM's own length (8 for BOS+6+EOS)
    is not what real_token_lengths returns.
    """
    residues = 6
    esm_tokenized_length = residues + 2
    assert real_token_lengths([residues + 1]) == [residues]
    assert real_token_lengths([residues + 1]) != [esm_tokenized_length]


if __name__ == "__main__":
    test_rare_residues_become_x()
    test_residues_are_whitespace_separated()
    test_eos_is_excluded_from_the_pooled_mean()
    test_pooled_length_never_reaches_zero()
    test_esm_and_prott5_pool_different_spans_on_purpose()
    print("OK: ProtT5 sequence prep and pooled-span arithmetic")
