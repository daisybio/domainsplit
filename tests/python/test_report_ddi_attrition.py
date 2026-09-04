#!/usr/bin/env python3
"""Checks for bin/report_ddi_attrition.py.

The point of the report is that a source's DDIs are accounted for end to end, so
what is asserted here is arithmetic rather than formatting:

* ``in_db == pruned + surviving`` for every token;
* a comma-joined ``source`` increments *every* one of its tokens, so token rows
  overlap and only the ``ALL_SOURCES`` row counts distinct DDIs;
* ``to_splitting`` is attached to the ``--splitting-source`` token and is ``NA``
  everywhere else -- external sources are never split, and a number there would
  read as leakage;
* a token that exists only in the pruned report (every one of its DDIs is gone)
  still gets a row, because "all of it went" is the answer the report exists to
  give.

Run directly (`python3 tests/python/test_report_ddi_attrition.py`) or via pytest.
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

REPORT = os.path.join(BIN, "report_ddi_attrition.py")

from domainsplit_schema import add_instance, make_db  # noqa: E402

ENV = dict(os.environ, PYTHONPATH=BIN + os.pathsep + os.environ.get("PYTHONPATH", ""))


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def write_csv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)


def build_master(path):
    """Surviving DDIs: two 3did, one PPIDM contributed under two tokens, one negatome.

    Every family also gets one instance, so the tier breakdown has something to
    count. The strata are deliberately mixed: PF00003 rests on a non-human
    reviewed protein and PF00005 on a human TrEMBL one, so the DDIs touching them
    must land above ``human_reviewed`` -- those are exactly the DDIs a narrower
    ``--instance_tier`` would have pruned.
    """
    make_db(path)
    conn = sqlite3.connect(path)
    conn.execute("PRAGMA foreign_keys=ON")
    pfams = ["PF00001", "PF00002", "PF00003", "PF00004", "PF00005", "PF00006"]
    conn.executemany("INSERT INTO domain(pfam_id) VALUES (?)", [(p,) for p in pfams])
    ids = {p: i for p, i in conn.execute("SELECT pfam_id, id FROM domain")}
    conn.executemany(
        "INSERT INTO domain_domain_interaction(domain_id_a, domain_id_b, negative, source) "
        "VALUES (?, ?, ?, ?)",
        [
            (ids["PF00001"], ids["PF00002"], 0, "3did"),
            (ids["PF00001"], ids["PF00003"], 0, "3did"),
            (ids["PF00004"], ids["PF00005"], 0, "PPIDM,PPIDM_Gold"),
            (ids["PF00004"], ids["PF00006"], 1, "negatome"),
        ],
    )
    # (family, protein, taxon, reviewed) -> stratum
    #   PF00001/2/4/6  human   reviewed    -> human_reviewed
    #   PF00003        mouse   reviewed    -> other_reviewed
    #   PF00005        human   unreviewed  -> human_unreviewed
    for pfam, uniprot, taxon, reviewed in [
        ("PF00001", "P00001", "9606", "reviewed"),
        ("PF00002", "P00002", "9606", "reviewed"),
        ("PF00003", "P00003", "10090", "reviewed"),
        ("PF00004", "P00004", "9606", "reviewed"),
        ("PF00005", "P00005", "9606", "unreviewed"),
        ("PF00006", "P00006", "9606", "reviewed"),
    ]:
        add_instance(conn, pfam, uniprot, 1, 10, taxon=taxon, reviewed=reviewed)
    conn.commit()
    conn.close()


def run(tmp):
    db = os.path.join(tmp, "domainsplit.sqlite3")
    build_master(db)

    pruned = os.path.join(tmp, "pruned_ddis.tsv")
    write_tsv(
        pruned,
        ["pfam_a", "pfam_b", "negative", "source", "reason"],
        [
            # 3did loses three; PPIDM loses one under both of its tokens;
            # single_domain_ppi loses everything it ever had.
            ("PF00007", "PF00008", 0, "3did", "no_instances:PF00007"),
            ("PF00007", "PF00009", 0, "3did", "no_instances:PF00007"),
            ("PF00010", "PF00011", 0, "3did", "no_instances:PF00010"),
            ("PF00012", "PF00013", 0, "PPIDM,PPIDM_Gold", "no_instances:PF00012"),
            ("PF00014", "PF00015", 0, "single_domain_ppi", "no_instances:PF00014"),
        ],
    )

    split_ddis = os.path.join(tmp, "split_ddis.csv")
    write_csv(split_ddis, ["protein1", "protein2"], [("PF00001", "PF00002")])

    external = os.path.join(tmp, "external.tsv")
    write_tsv(
        external,
        ["pfam_a", "pfam_b", "negative", "source"],
        [
            ("PF00004", "PF00005", 0, "PPIDM"),
            ("PF00004", "PF00005", 0, "PPIDM_Gold"),
            ("PF00012", "PF00013", 0, "PPIDM"),
            ("PF00012", "PF00013", 0, "PPIDM_Gold"),
            ("PF00014", "PF00015", 0, "single_domain_ppi"),
            ("PF00004", "PF00006", 1, "negatome"),
            ("PF00016", "PF00017", 1, "negatome"),
        ],
    )

    conflicts = os.path.join(tmp, "source_conflicts.tsv")
    write_tsv(
        conflicts,
        ["pfam_a", "pfam_b", "source", "negative", "action"],
        [("PF00016", "PF00017", "negatome", 1, "dropped_label_conflict")],
    )

    counts = os.path.join(tmp, "insert_3did_counts.tsv")
    write_tsv(counts, ["source", "offered", "inserted"], [("3did", 9, 5)])

    out = os.path.join(tmp, "ddi_source_attrition.tsv")
    tier_out = os.path.join(tmp, "ddi_tier_breakdown.tsv")
    subprocess.run(
        [
            sys.executable, REPORT,
            "--db", db,
            "--pruned-ddis", pruned,
            "--split-ddis", split_ddis,
            "--external-ddis", external,
            "--conflicts", conflicts,
            "--offered-counts", counts,
            "--splitting-source", "3did",
            "--out", out,
            "--tier-out", tier_out,
            "--versions", os.path.join(tmp, "versions.yml"),
            "--process-name", "TEST:REPORT_DDI_ATTRITION",
        ],
        check=True,
        env=ENV,
        capture_output=True,
        text=True,
    )

    with open(out) as fh:
        body = [line for line in fh if not line.startswith("#")]
    rows = {r["source"]: r for r in csv.DictReader(body, delimiter="\t")}

    with open(tier_out) as fh:
        tier_body = [line for line in fh if not line.startswith("#")]
    tiers = {
        (r["scope"], r["tier"]): int(r["count"])
        for r in csv.DictReader(tier_body, delimiter="\t")
    }
    return rows, tiers


def test_report():
    with tempfile.TemporaryDirectory() as tmp:
        rows, _tiers = run(tmp)

    # 3did: 9 offered, 5 inserted (dedup), 3 pruned, 2 surviving, 1 exported.
    assert rows["3did"]["offered"] == "9"
    assert rows["3did"]["in_db"] == "5"
    assert rows["3did"]["to_splitting"] == "1"
    assert rows["3did"]["pruned"] == "3"
    assert rows["3did"]["surviving"] == "2"

    # One DDI, two tokens: both are counted, so the two rows are identical and
    # neither may be added to the other.
    for token in ("PPIDM", "PPIDM_Gold"):
        assert rows[token]["offered"] == "2", token
        assert rows[token]["pruned"] == "1", token
        assert rows[token]["surviving"] == "1", token
        assert rows[token]["in_db"] == "2", token
        assert rows[token]["to_splitting"] == "NA", token

    # A source that lost everything still has a row.
    assert rows["single_domain_ppi"]["surviving"] == "0"
    assert rows["single_domain_ppi"]["pruned"] == "1"

    # The conflict column is not a subset of in_db: that pair never landed.
    assert rows["negatome"]["offered"] == "2"
    assert rows["negatome"]["in_db"] == "1"
    assert rows["negatome"]["conflict"] == "1"

    # in_db == pruned + surviving, for every token.
    for source, row in rows.items():
        if source == "ALL_SOURCES":
            continue
        assert int(row["in_db"]) == int(row["pruned"]) + int(row["surviving"]), source

    # ALL_SOURCES counts distinct DDI rows, so it is *less* than the token sums
    # wherever a pair carries several tokens.
    total = rows["ALL_SOURCES"]
    assert total["surviving"] == "4"
    assert total["pruned"] == "5"
    assert total["in_db"] == "9"
    assert total["to_splitting"] == "1"
    assert total["conflict"] == "1"
    assert sum(int(r["surviving"]) for s, r in rows.items() if s != "ALL_SOURCES") > 4


def test_tier_breakdown():
    """The second report: which stratum each surviving DDI actually rests on.

    A DDI is filed under the *worse* of its two families' best strata, which is
    the "would a narrower --instance_tier have kept this" question. With the
    fixture's instances:

      PF00001-PF00002  human_reviewed  x human_reviewed    -> human_reviewed
      PF00001-PF00003  human_reviewed  x other_reviewed    -> other_reviewed
      PF00004-PF00005  human_reviewed  x human_unreviewed  -> human_unreviewed
      PF00004-PF00006  human_reviewed  x human_reviewed    -> human_reviewed
    """
    with tempfile.TemporaryDirectory() as tmp:
        _rows, tiers = run(tmp)

    assert tiers[("ddis", "human_reviewed")] == 2
    assert tiers[("ddis", "other_reviewed")] == 1
    assert tiers[("ddis", "human_unreviewed")] == 1
    assert ("ddis", "other_unreviewed") not in tiers
    # Nothing may fall through the cracks: PRUNE_UNREPRESENTED_DDIS guarantees
    # every surviving DDI has instances on both sides.
    assert ("ddis", "no_instance") not in tiers
    assert sum(n for (scope, _t), n in tiers.items() if scope == "ddis") == 4

    # Families and proteins are counted under their own best stratum, one each.
    assert tiers[("families", "human_reviewed")] == 4
    assert tiers[("families", "other_reviewed")] == 1
    assert tiers[("families", "human_unreviewed")] == 1
    assert tiers[("proteins", "other_reviewed")] == 1
    assert tiers[("instances", "human_unreviewed")] == 1


if __name__ == "__main__":
    test_report()
    test_tier_breakdown()
    print("ok")
