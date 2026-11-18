#!/usr/bin/env python3
"""Local unit-check for bin/build_ppidm_sources.py on tiny synthetic inputs."""

import csv
import gzip
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BIN = os.path.join(REPO, "bin")
SCRIPT = os.path.join(BIN, "build_ppidm_sources.py")


def entry(id_, status, accs, taxid, pfams=(), strings=()):
    lines = [f"ID   {id_}   {status};         100 AA.",
             f"AC   {'; '.join(accs)};",
             f"OX   NCBI_TaxID={taxid} {{ECO:0000313}};"]
    lines += [f"DR   Pfam; {p}; Name; 1." for p in pfams]
    lines += [f"DR   STRING; {s}; -." for s in strings]
    # sequence block: a Pfam-looking line here must not be read
    lines += ["SQ   SEQUENCE   100 AA;", "     DR   Pfam; PF09999; X; 1.", "//"]
    return "\n".join(lines)


DAT = "\n".join([
    # O00001: secondary of P11111 only -> resolves; O00002: secondary of two entries -> does not
    entry("A_HUMAN", "Reviewed", ["P11111", "O00001"], 9606, ["PF00002", "PF00001", "PF00002"], ["9606.ENSP1"]),
    entry("B_HUMAN", "Reviewed", ["Q22222", "O00002"], 9606, ["PF00003"], ["9606.ENSP2"]),
    entry("C_HUMAN", "Reviewed", ["P44444", "O00002"], 9606, ["PF00004"], ["9606.ENSP4", "9606.ENSPAMB"]),
    entry("D_HUMAN", "Reviewed", ["P55555"], 9606, ["PF00005"], ["9606.ENSP5", "9606.ENSPAMB"]),
    entry("E_MOUSE", "Reviewed", ["M66666"], 10090, ["PF00006"], ["10090.ENSMUSP6"]),
    entry("F_HUMAN", "Unreviewed", ["A77777", "O00003"], 9606, ["PF00007"]),
    entry("G_HUMAN", "Reviewed", ["R88888"], 9606),  # no Pfam: outside the universe
]) + "\n"

EXTRACT_HEADER = "acc_a\tacc_b\tsource_db\tscore\tpmids"
BIOGRID = [
    "P11111\tQ22222\tbiogrid\t-\t100;200",
    "O00001\tQ22222\tbiogrid\t-\t200;300",  # secondary -> same pair, 3 distinct pubs
    "P44444\tP55555\tbiogrid\t-\t100",
    "P44444\tP55555\tbiogrid\t-\t100",      # same pub twice -> 1 distinct
    "O00002\tP11111\tbiogrid\t-\t1",        # shared secondary -> unmapped
    "A77777\tP11111\tbiogrid\t-\t1",        # unreviewed -> unmapped
    "O00003\tP11111\tbiogrid\t-\t1",        # secondary of unreviewed -> unmapped
    "O00001\tP11111\tbiogrid\t-\t1",        # self after mapping
]
INTACT = [
    "P11111\tQ22222\t0469\t0.3\t1",
    "P11111\tQ22222\t0469\t0.45\t1",  # intact max 0.45 > 0.4
    "P44444\tP55555\t0469\t0.4\t1",   # not strictly above
    "P44444\tP55555\t0471\t0.9\t1",   # mint
    "P11111\tP44444\t0465\t-\t1",     # dip, no score
    "P11111\tR88888\t0469\t0.9\t1",   # unmapped
]
HIPPIE = [
    "uniprot_accession_A\tuniprot_name_A\tentrez_A\tuniprot_accession_B\tuniprot_name_B\tentrez_B\tscore\tinfo",
    "P11111\tA\t1\tQ22222\tB\t2\t0.72\tx",
    "Q22222\tB\t2\tP11111\tA\t1\t0.8000\tx",  # reverse direction, higher score
    "P44444\tC\t4\tP55555\tD\t5\t0.72\tx",    # not strictly above
    "P11111-2\tA\t1\tP55555\tD\t5\t0.5\tx",   # isoform suffix
    "\t-\t9\tP55555\tD\t5\t0.9\tx",           # no accession -> unmapped
]
STRING_HEADER = "protein1 protein2 neighborhood fusion cooccurence coexpression experimental database textmining combined_score"
STRING_HUMAN = [
    "9606.ENSP1 9606.ENSP2 0 0 0 0 500 0 0 950",
    "9606.ENSP2 9606.ENSP1 0 0 0 0 500 0 0 950",
    "9606.ENSP4 9606.ENSP1 0 0 0 0 0 0 0 900",
    "9606.ENSP1 9606.ENSP4 0 0 0 0 0 0 0 900",
    "9606.ENSPAMB 9606.ENSP1 0 0 0 0 300 0 0 300",  # ambiguous id -> unmapped
    "9606.ENSP5 9606.ENSP1 0 0 0 0 300 0 0 300",    # experimental only
]
STRING_MOUSE = ["10090.ENSMUSP6 10090.ENSMUSPX 0 0 0 0 0 0 0 500"]
SIFTS = [
    "# 2026/11/08 - 12:37 | PDB: 45.26 | UniProt: 2026.04",
    "PDB\tCHAIN\tSP_PRIMARY\tRES_BEG\tRES_END\tPDB_BEG\tPDB_END\tSP_BEG\tSP_END",
    "1ABC\tA\tP11111\t1\t9\t1\t9\t1\t9",
    "1abc\tB\tQ22222\t1\t9\t1\t9\t1\t9",
    "1abc\tB\tP44444-2\t1\t9\t1\t9\t1\t9",
    "2xyz\tA\tP11111\t1\t9\t1\t9\t1\t9",
    "2xyz\tB\tM66666\t1\t9\t1\t9\t1\t9",   # mouse: cross-taxon
    "3aaa\tA\tP55555\t1\t9\t1\t9\t1\t9",
    "3aaa\tC\tR88888\t1\t9\t1\t9\t1\t9",   # no Pfam -> chain unmapped
    "4bbb\tA\tP44444\t1\t9\t1\t9\t1\t9",
    "4bbb\tB\tP44444\t1\t9\t1\t9\t1\t9",   # homodimer -> no pair
    "5nna\tNA\tP55555\t1\t9\t1\t9\t1\t9",  # chain literally "NA"
    "5nna\tA\tQ22222\t1\t9\t1\t9\t1\t9",
]
DID = [
    "#=ID\tPF00001\tPF00003\t (PF00001.1@Pfam PF00003.1@Pfam)",
    "#=3D\t1abc\tA:1-10\tB:1-10\t1.0\t1.0\t1:1",
    "#=3D\t1abc\tB:1-10\tA:1-10\t1.0\t1.0\t1:1",  # same chain pair reversed
    "#=3D\t1abc\tA:1-10\tA:20-30\t1.0\t1.0\t1:1",  # same chain
    "#=3D\t2xyz\tA:1-10\tB:1-10\t1.0\t1.0\t1:1",
    "#=3D\t3aaa\tA:1-10\tC:1-10\t1.0\t1.0\t1:1",
    "#=3D\t4bbb\tA:1-10\tB:1-10\t1.0\t1.0\t1:1",
    "#=3D\t5nna\tNA:1-10\tA:1-10\t1.0\t1.0\t1:1",
]


def write(path, lines, gz):
    with (gzip.open if gz else open)(path, "wt") as f:
        f.write("\n".join(lines) + "\n")
    return path


def read(path):
    with (gzip.open if path.endswith(".gz") else open)(path, "rt") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def pairs(outdir, name):
    return {(r["acc_a"], r["acc_b"]): int(r["hi_conf"]) for r in read(os.path.join(outdir, f"source_{name}.tsv.gz"))}


def run(tmp, hippie_gz):
    with gzip.open(os.path.join(tmp, "sprot.dat.gz"), "wt") as f:
        f.write(DAT)
    out = os.path.join(tmp, f"out_{hippie_gz}")
    cmd = [sys.executable, SCRIPT, "--sprot-dat", os.path.join(tmp, "sprot.dat.gz"),
           "--biogrid", write(os.path.join(tmp, "bg.tsv.gz"), [EXTRACT_HEADER] + BIOGRID, True),
           "--intact", write(os.path.join(tmp, "ia.tsv.gz"), [EXTRACT_HEADER] + INTACT, True),
           # misleading suffixes on purpose: format is sniffed, not read off the name
           "--hippie", write(os.path.join(tmp, "hippie.txt"), HIPPIE, hippie_gz),
           "--string-detailed", write(os.path.join(tmp, "9606.txt.gz"), [STRING_HEADER] + STRING_HUMAN, True),
           "--string-detailed", write(os.path.join(tmp, "10090.txt"), [STRING_HEADER] + STRING_MOUSE, False),
           "--sifts-chain-uniprot", write(os.path.join(tmp, "sifts.tsv.gz"), SIFTS, True),
           "--did-flat", write(os.path.join(tmp, "3did_flat.gz"), DID, True),
           "--hippie-min-score", "0.72", "--string-min-score", "900", "--intact-min-miscore", "0.4",
           "--biogrid-min-pubs", "2", "--outdir", out,
           "--versions", os.path.join(tmp, "v.yml"), "--process-name", "X"]
    r = subprocess.run(cmd, capture_output=True, text=True, env=dict(os.environ, PYTHONPATH=BIN))
    assert r.returncode == 0, r.stdout + r.stderr
    return out


def check(out):
    assert [(r["acc"], r["pfams"]) for r in read(os.path.join(out, "protein_pfam.tsv.gz"))] == [
        ("M66666", "PF00006"), ("P11111", "PF00001;PF00002"),
        ("P44444", "PF00004"), ("P55555", "PF00005"), ("Q22222", "PF00003")]

    assert pairs(out, "biogrid") == {("P11111", "Q22222"): 1, ("P44444", "P55555"): 0}
    assert pairs(out, "intact") == {("P11111", "Q22222"): 1, ("P44444", "P55555"): 0}
    assert pairs(out, "mint") == {("P44444", "P55555"): 1}
    assert pairs(out, "dip") == {("P11111", "P44444"): 0}
    assert pairs(out, "hippie") == {("P11111", "Q22222"): 1, ("P44444", "P55555"): 0, ("P11111", "P55555"): 0}
    assert pairs(out, "string_exp") == {("P11111", "Q22222"): 1, ("P11111", "P55555"): 0}
    assert pairs(out, "string_rest") == {("P11111", "Q22222"): 1, ("P11111", "P44444"): 0}
    assert pairs(out, "sifts") == {("P11111", "Q22222"): 1, ("P11111", "P44444"): 1, ("P55555", "Q22222"): 1}

    stats = {r["source"]: r for r in read(os.path.join(out, "ppidm_sources.tsv"))}
    assert list(stats) == ["intact", "mint", "dip", "biogrid", "hippie", "string_exp", "string_rest", "sifts"]
    assert stats["biogrid"] == {"source": "biogrid", "raw_rows": "8", "mapped_pairs": "2",
                                "hi_conf_pairs": "1", "unmapped_rows": "3"}
    assert stats["intact"]["raw_rows"] == "4" and stats["intact"]["unmapped_rows"] == "1"
    assert stats["hippie"]["raw_rows"] == "5" and stats["hippie"]["unmapped_rows"] == "1"
    assert (stats["string_exp"]["raw_rows"], stats["string_exp"]["unmapped_rows"]) == ("4", "1")
    assert (stats["string_rest"]["raw_rows"], stats["string_rest"]["unmapped_rows"]) == ("5", "1")
    assert (stats["sifts"]["raw_rows"], stats["sifts"]["unmapped_rows"]) == ("5", "1")


def test_build_ppidm_sources():
    with tempfile.TemporaryDirectory() as tmp:
        outs = [run(tmp, gz) for gz in (True, False)]
        for out in outs:
            check(out)
        assert read(os.path.join(outs[0], "source_hippie.tsv.gz")) == read(os.path.join(outs[1], "source_hippie.tsv.gz"))


if __name__ == "__main__":
    test_build_ppidm_sources()
    print("ok")
