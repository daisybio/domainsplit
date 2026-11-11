#!/usr/bin/env python3
"""fetch_afdb_structures.py's placement and mismatch check, on synthetic chains.

AFDB numbers residues by UniProt position, but the model was predicted from an
older sequence version. These pin that a residue lands at its own UniProt
position (a chain whose numbering starts mid-protein is placed there, not at
index 0), that a sequence that has since changed is measured as mismatches,
and that low-pLDDT residues are NaN rather than placed.

Run directly (`python3 tests/python/test_afdb_structures.py`) or via pytest.
"""

import gzip
import os
import sys
import tempfile

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "bin"))

from fetch_afdb_structures import (  # noqa: E402
    accessions_from_instance_fasta,
    backbone_resolved,
    place_on_uniprot,
    uniprot_sequences,
)

UNIPROT = "MKTAYIAKQRQISFVKSHFSRQ"  # L = 22


def chain(first_pos, seq):
    """A chain covering UniProt positions first_pos.. ; atom x-coordinate = position."""
    idx = list(range(first_pos, first_pos + len(seq)))
    atom37 = np.zeros((len(seq), 37, 3), dtype=np.float32)
    atom37[:, :, 0] = np.array(idx)[:, None]
    return idx, seq, atom37, np.full(len(seq), 90.0)


def test_numbering_offset_places_residues_at_their_uniprot_position():
    idx, seq, atom37, plddt = chain(5, UNIPROT[4:15])  # positions 5..15
    coords, placed, rate = place_on_uniprot(idx, seq, atom37, plddt, UNIPROT, 70)
    assert coords.shape == (22, 37, 3) and placed == 11 and rate == 0.0
    resolved = backbone_resolved(coords)
    assert resolved.tolist() == [False] * 4 + [True] * 11 + [False] * 7
    assert coords[4, 1, 0] == 5, "UniProt position 5 is array index 4"


def test_shifted_sequence_is_counted_as_mismatch():
    # The same residues numbered one too far: what a changed N-terminus does.
    idx, seq, atom37, plddt = chain(6, UNIPROT[4:15])
    coords, placed, rate = place_on_uniprot(idx, seq, atom37, plddt, UNIPROT, 70)
    assert rate > 0.5
    # Mismatching residues are not placed even when the protein would pass.
    assert backbone_resolved(coords).sum() == placed - round(rate * placed)


def test_out_of_range_positions_count_as_mismatch():
    idx, seq, atom37, plddt = chain(20, UNIPROT[19:] + "AA")  # runs past L
    _coords, placed, rate = place_on_uniprot(idx, seq, atom37, plddt, UNIPROT, 70)
    assert placed == 5 and rate == 2 / 5


def test_low_plddt_residues_are_nan():
    idx, seq, atom37, plddt = chain(1, UNIPROT)
    plddt[3:6] = 40.0
    coords, _placed, rate = place_on_uniprot(idx, seq, atom37, plddt, UNIPROT, 70)
    assert rate == 0.0
    assert (~backbone_resolved(coords)).nonzero()[0].tolist() == [3, 4, 5]


def test_input_parsers():
    with tempfile.TemporaryDirectory() as tmp:
        dom = os.path.join(tmp, "d.fasta")
        with open(dom, "w") as f:
            f.write(">PF00001_P12345_3_9\nAAA\n>PF00002_P12345_20_30\nCC\n>PF00003_Q9XYZ1_1_5\nG\n")
        assert accessions_from_instance_fasta(dom) == {"P12345", "Q9XYZ1"}
        uni = os.path.join(tmp, "u.fasta.gz")
        with gzip.open(uni, "wt") as f:
            f.write(">sp|P12345|A_HUMAN desc OS=x OX=9606\nMKT\nAYI\n>tr|O00000|B desc\nGG\n")
        assert uniprot_sequences(uni, {"P12345"}) == {"P12345": "MKTAYI"}


if __name__ == "__main__":
    test_numbering_offset_places_residues_at_their_uniprot_position()
    test_shifted_sequence_is_counted_as_mismatch()
    test_out_of_range_positions_count_as_mismatch()
    test_low_plddt_residues_are_nan()
    test_input_parsers()
    print("ok")
