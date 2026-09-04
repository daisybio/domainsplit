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
]


def test_parse_negatome():
    with tempfile.TemporaryDirectory() as tmp:
        negatome = os.path.join(tmp, "combined_pfam.txt")
        with open(negatome, "w") as fh:
            fh.write("\n".join(NEGATOME_LINES) + "\n")
        out = os.path.join(tmp, "negatome_ddis.tsv")

        env = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))
        subprocess.run(
            [sys.executable, PARSER, "--negatome", negatome, "--out", out,
             "--versions", os.path.join(tmp, "versions.yml"),
             "--process-name", "TEST:PARSE_NEGATOME"],
            check=True, env=env,
        )

        with open(out) as fh:
            rows = [
                (r["pfam_a"], r["pfam_b"], int(r["negative"]), r["source"])
                for r in csv.DictReader(fh, delimiter="\t")
            ]

        assert rows == [
            ("PF00001", "PF00002", 1, "negatome"),
            ("PF00003", "PF00004", 1, "negatome"),
        ], rows


if __name__ == "__main__":
    test_parse_negatome()
    print("OK: parse_negatome invariants hold")
