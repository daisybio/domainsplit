#!/usr/bin/env python3
"""Local unit-check for bin/parse_ppidm.py (no Nextflow, no cluster).

Runs the parser against a small synthetic ``predicted_ddi_ppi.tsv``, asserting:

  * domain tokens like ``10114/PF00069`` are parsed down to the Pfam accession;
  * every kept DDI is emitted twice -- under the umbrella source ``PPIDM`` and
    under ``PPIDM_<Class>`` -- so a consumer can select all of PPIDM or one
    confidence class;
  * a pair appearing under two classes keeps *both* class labels (they merge onto
    one DB row later) and still contributes only one ``PPIDM`` row;
  * unparseable tokens are skipped, and ``--classes`` filters which classes are
    emitted at all.

Run directly (`python3 tests/python/test_parse_ppidm.py`) or via pytest.
"""

import csv
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
PARSER = os.path.join(BIN, "parse_ppidm.py")

# Tokens carry a leading numeric id before the slash, as in real PPIDM output.
PPIDM_ROWS = [
    "domain_1\tdomain_2\tclass",          # header (skipped)
    "10/PF00001\t20/PF00002\tGold",       # kept -> PPIDM + PPIDM_Gold
    "30/PF00003\t40/PF00004\tSilver",     # kept -> PPIDM + PPIDM_Silver
    "50/PF00005\t60/PF00006\tBronze",     # kept -> PPIDM + PPIDM_Bronze
    "10/PF00001\t20/PF00002\tSilver",     # same pair, second class
    "junk\tnonsense\tGold",               # unparseable -> skipped
]


def run_parser(ppidm, classes, out, tmp):
    env = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))
    subprocess.run(
        [sys.executable, PARSER, "--ppidm", ppidm, "--classes", classes,
         "--out", out,
         "--versions", os.path.join(tmp, "versions.yml"),
         "--process-name", "TEST:PARSE_PPIDM"],
        check=True, env=env,
    )


def read_tsv(path):
    with open(path) as fh:
        return [
            (r["pfam_a"], r["pfam_b"], int(r["negative"]), r["source"])
            for r in csv.DictReader(fh, delimiter="\t")
        ]


def write_ppidm(path, rows):
    with open(path, "w") as fh:
        fh.write("\n".join(rows) + "\n")


def test_parse_ppidm_all_classes():
    with tempfile.TemporaryDirectory() as tmp:
        ppidm = os.path.join(tmp, "predicted_ddi_ppi.tsv")
        out = os.path.join(tmp, "ppidm_ddis.tsv")
        write_ppidm(ppidm, PPIDM_ROWS)
        run_parser(ppidm, "Bronze,Silver,Gold", out, tmp)

        rows = read_tsv(out)
        assert all(neg == 0 for _, _, neg, _ in rows), "PPIDM rows must be positives"

        by_pair = {}
        for a, b, _, source in rows:
            by_pair.setdefault((a, b), []).append(source)

        # The doubled pair keeps both class labels but only one umbrella row.
        assert sorted(by_pair[("PF00001", "PF00002")]) == [
            "PPIDM", "PPIDM_Gold", "PPIDM_Silver"
        ], by_pair[("PF00001", "PF00002")]
        assert sorted(by_pair[("PF00003", "PF00004")]) == ["PPIDM", "PPIDM_Silver"]
        assert sorted(by_pair[("PF00005", "PF00006")]) == ["PPIDM", "PPIDM_Bronze"]
        assert len(by_pair) == 3, f"unparseable row leaked: {sorted(by_pair)}"

        # Gold is written before Silver, so the merged source list downstream is
        # ordered by confidence.
        sources = [s for _, _, _, s in rows]
        assert sources.index("PPIDM_Gold") < sources.index("PPIDM_Silver")


def test_parse_ppidm_class_filter():
    """--classes restricts which classes are emitted at all."""
    with tempfile.TemporaryDirectory() as tmp:
        ppidm = os.path.join(tmp, "predicted_ddi_ppi.tsv")
        out = os.path.join(tmp, "ppidm_ddis.tsv")
        write_ppidm(ppidm, PPIDM_ROWS)
        run_parser(ppidm, "Gold", out, tmp)

        rows = read_tsv(out)
        assert {s for _, _, _, s in rows} == {"PPIDM", "PPIDM_Gold"}
        assert {(a, b) for a, b, _, _ in rows} == {("PF00001", "PF00002")}


if __name__ == "__main__":
    test_parse_ppidm_all_classes()
    test_parse_ppidm_class_filter()
    print("OK: parse_ppidm class handling holds")
