#!/usr/bin/env python3
"""Checks for bin/export_split_inputs.py.

The two files this script writes sit on opposite sides of the Pfam fetch, and
each one enforces a rule the rest of the pipeline depends on:

``--mode families``
    the fetch must cover 3did *and* the external sources, because
    BUILD_EXTERNAL_TEST needs instances for external families too. Asserted:
    the union, not just what the DB holds; deduplicated; accession-sorted.

``--mode ddis``
    only 3did positives may enter the split population, and only families that
    actually resolved to an instance. Asserted: external and negative rows never
    appear; a DDI touching an instance-less family is dropped; output is the
    ``protein1,protein2`` header ppi-splitting's read_ppis() expects.

Run directly (`python3 tests/python/test_export_split_inputs.py`) or via pytest.
"""

import csv
import os
import sqlite3
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

EXPORT = os.path.join(BIN, "export_split_inputs.py")

from ddi_db_utils import ensure_domains, insert_ddis, merge_external_ddis  # noqa: E402
from external_ddi_tsv import write_external_tsv  # noqa: E402
from domainsplit_schema import add_instance, make_db  # noqa: E402

ENV = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))

INSTANCE_COLUMNS = ["instance_id", "family", "clan", "protein_id", "start", "end", "taxon_id", "source_db"]


def build_db(tmp):
    """3did over PF00002/PF00010/PF00001, plus external pairs reaching PF00300."""
    path = os.path.join(tmp, "domainsplit.sqlite3")
    make_db(path)
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA foreign_keys=ON")
    ensure_domains(conn, {"PF00001", "PF00002", "PF00010", "PF00300"})
    insert_ddis(
        conn,
        [("PF00002", "PF00010"), ("PF00001", "PF00002"), ("PF00010", "PF00010")],
        negative=False,
        source="3did",
    )
    # An external positive and an external negative. Neither may reach the
    # exported split population, whatever its label.
    merge_external_ddis(
        conn,
        [
            ("PF00001", "PF00300", 0, "PPIDM"),
            ("PF00002", "PF00300", 1, "negatome"),
        ],
    )
    conn.commit()
    conn.close()


def write_instances(path, families):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow(INSTANCE_COLUMNS)
        for i, family in enumerate(families):
            writer.writerow([f"{family}_P{i:05d}_1_50", family, "", f"P{i:05d}", 1, 50, "9606", "pfam"])


def run_export(tmp, mode, **kwargs):
    out = os.path.join(tmp, "out.txt")
    cmd = [
        sys.executable,
        EXPORT,
        "--mode",
        mode,
        "--db",
        os.path.join(tmp, "domainsplit.sqlite3"),
        "--out",
        out,
        "--versions",
        os.path.join(tmp, "versions.yml"),
        "--process-name",
        "TEST:EXPORT",
    ]
    for flag, value in kwargs.items():
        cmd.append(f"--{flag.replace('_', '-')}")
        cmd.extend(value if isinstance(value, list) else [value])
    result = subprocess.run(cmd, env=ENV, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return out, result.stdout


def test_families_is_the_union_sorted_and_deduplicated():
    with tempfile.TemporaryDirectory() as tmp:
        build_db(tmp)
        # PF00300 is already in the DB via the external merge; PF00500 is only in
        # the TSV. Both must appear exactly once.
        tsv = os.path.join(tmp, "extra.tsv")
        write_external_tsv(tsv, [("PF00500", "PF00300", 0, "single_domain_ppi")])

        out, _log = run_export(tmp, "families", external_ddis=[tsv])
        families = open(out).read().split()

        assert families == ["PF00001", "PF00002", "PF00010", "PF00300", "PF00500"]


def test_ddis_are_3did_positives_with_instances_only():
    with tempfile.TemporaryDirectory() as tmp:
        build_db(tmp)
        instances = os.path.join(tmp, "instances.tsv")
        # PF00001 resolved to no instance -- the instance_tier stranded it -- so the
        # 3did DDI touching it must not be offered for splitting.
        write_instances(instances, ["PF00002", "PF00010", "PF00300"])

        out, log = run_export(tmp, "ddis", instances=instances)
        with open(out) as fh:
            rows = list(csv.reader(fh))

        assert rows[0] == ["protein1", "protein2"]
        pairs = {tuple(row) for row in rows[1:]}
        # PF00002/PF00010 and the PF00010 self-pair survive; PF00001/PF00002 does
        # not, and no external pair appears at all.
        assert pairs == {("PF00002", "PF00010"), ("PF00010", "PF00010")}
        assert "PF00300" not in {p for pair in pairs for p in pair}
        assert "ddis_dropped_no_instances = 1" in log


def test_ddis_rejects_a_missing_instances_argument():
    with tempfile.TemporaryDirectory() as tmp:
        build_db(tmp)
        result = subprocess.run(
            [
                sys.executable,
                EXPORT,
                "--mode",
                "ddis",
                "--db",
                os.path.join(tmp, "domainsplit.sqlite3"),
                "--out",
                os.path.join(tmp, "out.csv"),
                "--versions",
                os.path.join(tmp, "versions.yml"),
                "--process-name",
                "TEST:EXPORT",
            ],
            env=ENV,
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "--instances" in result.stderr


if __name__ == "__main__":
    test_families_is_the_union_sorted_and_deduplicated()
    test_ddis_are_3did_positives_with_instances_only()
    test_ddis_rejects_a_missing_instances_argument()
    print("ok")
