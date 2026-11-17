#!/usr/bin/env python3
"""Local unit-check for bin/extract_ppi_mitab.py on synthetic BioGRID/IntAct MITAB."""

import csv
import gzip
import os
import subprocess
import sys
import tempfile
import zipfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
SCRIPT = os.path.join(BIN, "extract_ppi_mitab.py")

PHYS = 'psi-mi:"MI:0915"(physical association)'
GENETIC = 'psi-mi:"MI:2373"(negative genetic interaction)'
Y2H = 'psi-mi:"MI:0018"(two hybrid)'


def biogrid_row(alt_a, alt_b, tax_a="taxid:9606", tax_b="taxid:9606", itype=PHYS,
                pubs="pubmed:200|pubmed:100"):
    return "\t".join(["biogrid:1", "biogrid:2", alt_a, alt_b, "-", "-", Y2H, "Smith (2000)", pubs,
                      tax_a, tax_b, itype, 'psi-mi:"MI:0463"(biogrid)', "biogrid:9", "-"])


def intact_row(a, b, srcdb="MI:0469", score="intact-miscore:0.56", negative="false",
               tax_a="taxid:9606(human)|taxid:9606(Homo sapiens)", tax_b=None):
    c = ["-"] * 42
    c[0], c[1] = a, b
    c[6] = Y2H
    c[8] = "pubmed:300|imex:IM-1"
    c[9], c[10] = tax_a, tax_b or tax_a
    c[11] = PHYS
    c[12] = f'psi-mi:"{srcdb}"(x)'
    c[14] = f"author score:high|{score}"
    c[35] = negative
    return "\t".join(c)


BIOGRID = [
    "#ID Interactor A\t...",
    biogrid_row("entrez gene/locuslink:1|uniprot/swiss-prot:Q22222", "uniprot/swiss-prot:P11111"),  # kept, swapped
    biogrid_row("uniprot/swiss-prot:P11111", "uniprot/swiss-prot:Q22222", itype=GENETIC),           # genetic
    biogrid_row("uniprot/swiss-prot:P11111", "refseq:NP_1"),                                         # no swiss-prot
    biogrid_row("uniprot/swiss-prot:P11111|uniprot/swiss-prot:P33333", "uniprot/swiss-prot:Q22222"),  # several
    biogrid_row("uniprot/swiss-prot:P11111", "uniprot/swiss-prot:Q22222", tax_b="taxid:10090"),     # cross-taxon
    biogrid_row("uniprot/swiss-prot:P11111-2", "uniprot/swiss-prot:P11111"),                        # self after strip
    biogrid_row("uniprot/swiss-prot:P44444", "uniprot/swiss-prot:P55555", pubs="-"),                # kept, no pmid
]

INTACT = [
    "#ID(s) interactor A\t...",
    intact_row("uniprotkb:P11111-PRO_0000000001", "uniprotkb:Q22222-2",
               score="intact-miscore:0.3|intact-miscore:0.45"),                # kept, max score
    intact_row("uniprotkb:P11111", "chebi:\"CHEBI:1\""),                        # non-uniprot
    intact_row("uniprotkb:P11111", "uniprotkb:Q22222", negative="true"),        # negative
    intact_row("uniprotkb:P11111", "uniprotkb:Q22222", tax_a="taxid:-1(in vitro)"),  # missing taxon
    intact_row("uniprotkb:P11111", "uniprotkb:Q22222", tax_b="taxid:10090(mouse)"),  # cross-taxon
    intact_row("uniprotkb:P44444", "uniprotkb:P55555", srcdb="MI:0471"),        # MINT
    intact_row("uniprotkb:P44444", "uniprotkb:P55555", srcdb="MI:0465", score="-"),  # DIP, no score
]


def run(inp, fmt, tmp):
    out = os.path.join(tmp, f"out_{fmt}.tsv.gz")
    env = dict(os.environ, PYTHONPATH=BIN)
    r = subprocess.run([sys.executable, SCRIPT, "--input", inp, "--format", fmt, "--out", out,
                        "--versions", os.path.join(tmp, "v.yml"), "--process-name", "X"],
                       capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stdout + r.stderr
    with gzip.open(out, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t")), r.stdout


def write_inputs(tmp, name, lines, member):
    text = "\n".join(lines) + "\n"
    plain = os.path.join(tmp, f"{name}.txt")
    with open(plain, "w") as f:
        f.write(text)
    gz = os.path.join(tmp, f"{name}.mitab.gz")
    with gzip.open(gz, "wt") as f:
        f.write(text)
    zp = os.path.join(tmp, f"{name}.zip")
    with zipfile.ZipFile(zp, "w") as z:
        z.writestr(member, text)
        if name == "intact":  # must never be read
            z.writestr("intact_negative.txt", intact_row("uniprotkb:P99999", "uniprotkb:Q99999") + "\n")
    return [plain, gz, zp]


def test_biogrid():
    with tempfile.TemporaryDirectory() as tmp:
        for inp in write_inputs(tmp, "biogrid", BIOGRID, "BIOGRID-ALL-5.0.1.mitab.txt"):
            rows, log = run(inp, "biogrid", tmp)
            assert [(r["acc_a"], r["acc_b"]) for r in rows] == [("P11111", "Q22222"), ("P44444", "P55555")]
            assert rows[0] == {"acc_a": "P11111", "acc_b": "Q22222", "taxid": "9606", "method": "0018",
                               "source_db": "biogrid", "score": "-", "pmids": "100;200"}
            assert rows[1]["pmids"] == ""
            for reason in ("non_physical_type", "no_swissprot_id", "several_swissprot_ids",
                           "cross_taxon", "self_pair"):
                assert f"'{reason}': 1" in log, (reason, log)


def test_intact():
    with tempfile.TemporaryDirectory() as tmp:
        for inp in write_inputs(tmp, "intact", INTACT, "intact.txt"):
            rows, log = run(inp, "intact", tmp)
            got = [(r["acc_a"], r["acc_b"], r["source_db"], r["score"]) for r in rows]
            assert got == [("P11111", "Q22222", "0469", "0.45"), ("P44444", "P55555", "0471", "0.56"),
                           ("P44444", "P55555", "0465", "-")], got
            assert rows[0]["pmids"] == "300" and rows[0]["taxid"] == "9606"
            for reason in ("non_uniprot_interactor", "negative", "taxid_missing_or_negative", "cross_taxon"):
                assert f"'{reason}': 1" in log, (reason, log)


if __name__ == "__main__":
    test_biogrid()
    test_intact()
    print("ok")
