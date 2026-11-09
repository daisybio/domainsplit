#!/usr/bin/env python3
"""Local unit-check for bin/parse_negatome.py (no Nextflow, no cluster).

Asserts that Negatome's whitespace-separated Pfam pairs become normalized rows
tagged ``negative=1, source='negatome'``, in canonical pair order, with short or
blank lines skipped and swapped duplicates collapsed.

Run directly (`python3 tests/python/test_parse_negatome.py`) or via pytest.
"""

import csv
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
PARSER = os.path.join(BIN, "parse_negatome.py")

NEGATOME_LINES = [
    "PF00002 PF00001",       # kept, reordered to (PF00001, PF00002)
    "PF00001\tPF00002",      # same pair again -> collapsed
    "PF00003 PF00004 extra",  # trailing columns ignored
    "PF00005",               # too short -> skipped
    "",                      # blank -> skipped
    "PF00006.12 PF00007.3",  # version suffixes stripped
    "<!DOCTYPE html>",       # not two Pfam accessions -> skipped
]


def run_parser(tmp, lines, check=True):
    negatome = os.path.join(tmp, "combined_pfam.txt")
    with open(negatome, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    out = os.path.join(tmp, "negatome_ddis.tsv")
    env = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))
    proc = subprocess.run(
        [sys.executable, PARSER, "--negatome", negatome, "--out", out,
         "--versions", os.path.join(tmp, "versions.yml"),
         "--process-name", "TEST:PARSE_NEGATOME"],
        check=check, env=env, capture_output=True, text=True,
    )
    return proc, out


def test_parse_negatome():
    with tempfile.TemporaryDirectory() as tmp:
        proc, out = run_parser(tmp, NEGATOME_LINES)
        assert "skipped 2 lines" in proc.stdout, proc.stdout

        with open(out) as fh:
            rows = [
                (r["pfam_a"], r["pfam_b"], int(r["negative"]), r["source"])
                for r in csv.DictReader(fh, delimiter="\t")
            ]

        assert rows == [
            ("PF00001", "PF00002", 1, "negatome"),
            ("PF00003", "PF00004", 1, "negatome"),
            ("PF00006", "PF00007", 1, "negatome"),
        ], rows


def test_an_html_error_page_is_fatal():
    """An error page served with HTTP 200 must not become a green, Negatome-less run."""
    page = ["<!DOCTYPE html>", "<html><head><title>404 Not Found</title></head>",
            "<h1>Not Found</h1>", "</html>"]
    with tempfile.TemporaryDirectory() as tmp:
        proc, _out = run_parser(tmp, page, check=False)
    assert proc.returncode != 0
    assert "holds no Pfam pair" in proc.stderr, proc.stderr


if __name__ == "__main__":
    test_parse_negatome()
    test_an_html_error_page_is_fatal()
    print("OK: parse_negatome invariants hold")
