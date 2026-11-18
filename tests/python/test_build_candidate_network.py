#!/usr/bin/env python3
"""Pin the rules in bin/build_candidate_network.py that define the ilp_candidates pool.

  * offline (``--mapping-in``): only 3did families, existing DDIs dropped,
    ``n_tested >= --min-n-tested`` kept (boundary inclusive);
  * UniProt lookup (mocked): a primary-name hit outranks an alias hit, every gene
    entry of a record counts, and the ``Link: rel="next"`` page is followed.

Run directly (`python3 tests/python/test_build_candidate_network.py`) or via pytest.
"""

import json
import os
import sqlite3
import subprocess
import sys
import tempfile

import pyarrow as pa
import pyarrow.parquet as pq

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
sys.path.insert(0, BIN)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import build_candidate_network as bcn  # noqa: E402
from ddi_db_utils import ensure_domains, insert_ddis  # noqa: E402
from domainsplit_schema import make_db  # noqa: E402


def test_offline_candidate_network():
    with tempfile.TemporaryDirectory() as tmp:
        db = os.path.join(tmp, "db.sqlite3")
        make_db(db)
        conn = sqlite3.connect(db)
        ensure_domains(conn, ["PF00001", "PF00002", "PF00003"])
        insert_ddis(conn, [("PF00001", "PF00002"), ("PF00003", "PF00003")], negative=False, source="3did")
        conn.commit()
        conn.close()

        mapping = os.path.join(tmp, "mapping.json")
        with open(mapping, "w") as fh:
            json.dump({
                "gene_to_uniprot": {"G1": "P1", "G2": "P2", "G3": "P3"},
                # PF00009 is not a 3did family
                "uniprot_to_pfams": {"P1": ["PF00001", "PF00009"], "P2": ["PF00002"], "P3": ["PF00003"]},
            }, fh)

        parquet = os.path.join(tmp, "ppi.parquet")
        pq.write_table(pa.table({
            "gene_name_bait": ["G1", "G1", "G2", "G2", "G1"],
            "gene_name_prey": ["G2", "G3", "G3", "G2", "G5"],
            "n_tested": [5, 5, 3, 2, 9],
        }), parquet)

        out = os.path.join(tmp, "candidate_network.csv")
        subprocess.run(
            [sys.executable, os.path.join(BIN, "build_candidate_network.py"),
             "--db", db, "--parquet", parquet, "--mapping-in", mapping,
             "--mapping-out", os.path.join(tmp, "out.json"),
             "--network-out", out, "--min-n-tested", "3"],
            check=True, capture_output=True,
        )
        with open(out) as fh:
            # G1-G2: existing 3did DDI; G2-G2: n_tested 2 < 3; G5: unmapped
            assert fh.read() == "protein1,protein2\nPF00001,PF00003\nPF00002,PF00003\n"


class FakeResponse:
    def __init__(self, results, next_url=None):
        self._results = results
        self.links = {"next": {"url": next_url}} if next_url else {}

    def raise_for_status(self):
        pass

    def json(self):
        return {"results": self._results}


def entry(acc, genes, pfams=()):
    return {
        "primaryAccession": acc,
        "genes": genes,
        "uniProtKBCrossReferences": [{"database": "Pfam", "id": p} for p in pfams],
    }


def test_alias_and_paging(monkeypatch):
    page1 = [
        entry("P00001", [{"geneName": {"value": "ABC1"}}], ["PF00001"]),
        entry("P00002", [{"geneName": {"value": "XYZ9"}, "synonyms": [{"value": "ABC1"}]}], ["PF00002"]),
        entry("P00004", [{"geneName": {"value": "Q1"}, "synonyms": [{"value": "GHI3"}]}]),
    ]
    page2 = [
        entry("P00003", [{"geneName": {"value": "OTHER"}}, {"geneName": {"value": "DEF2"}}], ["PF00003"]),
        entry("P00005", [{"geneName": {"value": "Q2"}, "synonyms": [{"value": "GHI3"}]}]),
    ]
    calls = []

    def fake_get(url, headers=None, params=None, timeout=None):
        calls.append((url, params))
        return FakeResponse(page1, "https://next.page/2") if len(calls) == 1 else FakeResponse(page2)

    monkeypatch.setattr(bcn.requests, "get", fake_get)
    monkeypatch.setattr(bcn.time, "sleep", lambda s: None)

    g2u, u2p = bcn.fetch_gene_mappings({"ABC1", "DEF2", "GHI3"}, batch_size=50)

    assert calls[1] == ("https://next.page/2", None)
    # ABC1: primary hit beats alias; DEF2: second gene entry; GHI3: two alias hits -> ambiguous
    assert g2u == {"ABC1": "P00001", "DEF2": "P00003"}
    assert u2p["P00003"] == {"PF00003"}


if __name__ == "__main__":
    import pytest
    sys.exit(pytest.main([__file__, "-q"]))
