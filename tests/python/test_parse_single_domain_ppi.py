#!/usr/bin/env python3
"""Checks for bin/parse_single_domain_ppi.py's HIPPIE layout detection.

HIPPIE is read by UniProt accession only, so a header without both accession
columns must be rejected -- a positional fallback reads entrez ids (or the score
column itself) as accessions, resolves nothing and reports zero DDIs with exit 0.

Run directly (`python3 tests/python/test_parse_single_domain_ppi.py`) or via pytest.
"""

import gzip
import json
import os
import subprocess
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO, "bin"))

from parse_single_domain_ppi import detect_layout  # noqa: E402

CURRENT = "\t".join([
    "uniprot_accession_A", "uniprot_name_A", "entrez_A",
    "uniprot_accession_B", "uniprot_name_B", "entrez_B", "score", "info",
]) + "\n"


def test_current_layout_is_located_by_name():
    assert detect_layout(CURRENT) == (0, 3, 6)


def test_header_without_accession_columns_is_fatal():
    """Names `score`, but no accession column: must fail, not fall back to offsets."""
    header = "\t".join(["name_A", "entrez_A", "name_B", "entrez_B", "score", "info"]) + "\n"
    with pytest.raises(SystemExit) as excinfo:
        detect_layout(header)
    assert "uniprot_accession_A / uniprot_accession_B" in str(excinfo.value)


def test_header_with_one_accession_column_is_fatal():
    header = "\t".join(["uniprot_accession_A", "entrez_A", "entrez_B", "score"]) + "\n"
    with pytest.raises(SystemExit) as excinfo:
        detect_layout(header)
    assert "uniprot_accession_B" in str(excinfo.value)


def test_header_without_score_is_fatal():
    with pytest.raises(SystemExit):
        detect_layout("AL1A1_HUMAN\t216\tAL1A1_HUMAN\t216\t0.76\tinfo\n")


@pytest.mark.parametrize("gz", [False, True])
def test_hippie_plain_and_gzip_read_alike(tmp_path, gz):
    """Format is sniffed from the magic bytes, so a gzip named `.txt` still reads."""
    hippie = tmp_path / "HIPPIE-current.txt"
    text = CURRENT + "P00001\tA_HUMAN\t1\tP00002\tB_HUMAN\t2\t0.9\tinfo\n"
    if gz:
        with gzip.open(hippie, "wt") as fh:
            fh.write(text)
    else:
        hippie.write_text(text)
    smap = tmp_path / "map.json"
    smap.write_text(json.dumps({"accession_to_pfams": {"P00001": ["PF00001"], "P00002": ["PF00002"]}}))
    out = tmp_path / "out.tsv"
    bin_dir = os.path.join(REPO, "bin")
    subprocess.run(
        [sys.executable, os.path.join(bin_dir, "parse_single_domain_ppi.py"),
         "--hippie", str(hippie), "--swissprot-map", str(smap), "--min-score", "0.5",
         "--out", str(out), "--versions", str(tmp_path / "v.yml"), "--process-name", "T"],
        check=True, env=dict(os.environ, PYTHONPATH=bin_dir),
    )
    assert out.read_text().splitlines()[1:] == ["PF00001\tPF00002\t0\tsingle_domain_ppi"]


if __name__ == "__main__":
    test_current_layout_is_located_by_name()
    test_header_without_accession_columns_is_fatal()
    test_header_with_one_accession_column_is_fatal()
    test_header_without_score_is_fatal()
    print("OK: HIPPIE layout detection")
