#!/usr/bin/env python3
"""Checks for bin/parse_swissprot_dat.py (no Nextflow, no cluster).

This parser replaced two ``rest.uniprot.org`` stream queries, which are throttled
to ~4 KB/s per connection and drop before an all-reviewed GO query finishes. What
matters is that the three files it carves out of the flat file are still the exact
shapes their consumers parse, because those consumers were left untouched:

  * ``insert_protein_go_terms.py`` reads ``Entry``/``Gene Ontology IDs`` with
    pandas and splits the GO column on ``"; "``;
  * ``build_swissprot_pfam_map.py`` reads four positional columns and splits the
    Pfam column on ``";"``, the gene column on whitespace;
  * ``insert_protein_sequences.py`` reads the FASTA with ``Bio.SeqIO`` and takes
    the accession from ``record.id.split("|")[1]``, so the ``sp|ACC|NAME`` header
    is load-bearing.

Plus the two things that make the flat file usable at all: only the *primary*
accession is emitted (an ``AC`` line lists secondary ones too, which the REST TSV
omitted), and ``--taxon-ids`` drops entries by ``OX`` so the human-only pipeline
does not parse and then discard 96% of SwissProt.

Run directly (`python3 tests/python/test_parse_swissprot_dat.py`) or via pytest.
"""

import gzip
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
PARSER = os.path.join(REPO, "bin", "parse_swissprot_dat.py")

# One human entry with everything on it, one mouse entry the taxon filter must
# drop, one human entry with neither GO nor Pfam.
DAT = """\
ID   ARVC_HUMAN              Reviewed;         30 AA.
AC   O00192; A6NDS4; Q5T9Y1;
DE   RecName: Full=Armadillo repeat protein {ECO:0000305};
DE            Short=ARVCF;
GN   Name=ARVCF; Synonyms=ARVC, XYZ1 {ECO:0000312};
OS   Homo sapiens (Human).
OX   NCBI_TaxID=9606 {ECO:0000313|EMBL:AAA};
DR   GO; GO:0005856; C:cytoskeleton; IEA:UniProtKB-KW.
DR   GO; GO:0005515; F:protein binding; IPI:IntAct.
DR   Pfam; PF00514; Arm; 6.
DR   Pfam; PF00071; Ras; 1.
DR   STRING; 9606.ENSP00000263045; -.
DR   PROSITE; PS50176; ARM_REPEAT; 1.
SQ   SEQUENCE   30 AA;  3466 MW;  A1B2C3D4E5F6 CRC64;
     MSSRSVSRSR ARSRSPRRSR SRSRSPAYKR
//
ID   MOUSY_MOUSE             Reviewed;         12 AA.
AC   P99999;
DE   RecName: Full=Mouse only protein;
GN   Name=Mus1;
OS   Mus musculus (Mouse).
OX   NCBI_TaxID=10090;
DR   GO; GO:0000001; C:nope; IEA:X.
DR   Pfam; PF99999; Nope; 1.
DR   STRING; 10090.ENSMUSP00000000001; -.
SQ   SEQUENCE   12 AA;  1000 MW;  FFFF CRC64;
     MSSRSVSRSR AR
//
ID   NOGO_HUMAN              Reviewed;         10 AA.
AC   Q11111;
DE   RecName: Full=No GO no Pfam protein;
GN   OrderedLocusNames=b0001; ORFNames=CG42;
OS   Homo sapiens (Human).
OX   NCBI_TaxID=9606;
SQ   SEQUENCE   10 AA;  900 MW;  EEEE CRC64;
     MSSRSVSRSR
//
"""


def run(tmp, taxon_ids):
    dat = os.path.join(tmp, "sprot.dat.gz")
    with gzip.open(dat, "wt") as fh:
        fh.write(DAT)
    out = {name: os.path.join(tmp, name) for name in
           ("go.tsv.gz", "pfam.tsv.gz", "seq.fasta.gz", "string.tsv.gz")}
    proc = subprocess.run(
        [sys.executable, PARSER, "--dat", dat,
         "--go-terms", out["go.tsv.gz"],
         "--pfam-map", out["pfam.tsv.gz"],
         "--sequences", out["seq.fasta.gz"],
         "--string-map", out["string.tsv.gz"],
         "--taxon-ids", taxon_ids,
         "--versions", os.path.join(tmp, "versions.yml"),
         "--process-name", "TEST"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    return {name: gzip.open(path, "rt").read() for name, path in out.items()}


def test_string_map_replaces_the_idmapping_download():
    """The STRING map must be byte-compatible with what insert_ppi.py already reads.

    That script's loader splits each line on tabs into (uniprot, id_type, symbol)
    and keeps the rows whose type is "STRING", so the map is a drop-in for the
    per-organism <ORG>_<taxid>_idmapping.dat.gz it replaces -- and unlike that file
    it follows --taxon-ids, so it covers exactly the species the run's protein
    universe does.
    """
    with tempfile.TemporaryDirectory() as tmp:
        human = run(tmp, "9606")["string.tsv.gz"].splitlines()
        assert human == ["O00192\tSTRING\t9606.ENSP00000263045"], human

    with tempfile.TemporaryDirectory() as tmp:
        every = run(tmp, "")["string.tsv.gz"].splitlines()
        assert every == [
            "O00192\tSTRING\t9606.ENSP00000263045",
            "P99999\tSTRING\t10090.ENSMUSP00000000001",
        ], every

    # insert_ppi.py's own loader, run over the file it will actually be handed.
    sys.path.insert(0, os.path.dirname(PARSER))
    from insert_ppi import load_uniprot_id_mapping  # noqa: E402

    with tempfile.TemporaryDirectory() as tmp:
        run(tmp, "")
        mapping = load_uniprot_id_mapping(os.path.join(tmp, "string.tsv.gz"), "STRING")
    assert mapping == {
        "9606.ENSP00000263045": "O00192",
        "10090.ENSMUSP00000000001": "P99999",
    }, mapping


def test_human_filter_and_output_shapes():
    with tempfile.TemporaryDirectory() as tmp:
        got = run(tmp, "9606")

        go = got["go.tsv.gz"].splitlines()
        assert go[0] == "Entry\tGene Ontology IDs"
        # Primary accession only, GO ids joined the way pandas' consumer splits
        # them, evidence codes stripped, IEA kept.
        assert go[1] == "O00192\tGO:0005856; GO:0005515"
        # No GO line at all means no row -- not a row with an empty column.
        assert len(go) == 2, go

        pfam = got["pfam.tsv.gz"].splitlines()
        assert pfam[0] == "Entry\tEntry Name\tGene Names\tPfam"
        # Evidence braces stripped, comma-separated synonyms split, order kept.
        assert pfam[1] == "O00192\tARVC_HUMAN\tARVCF ARVC XYZ1\tPF00514;PF00071;"
        # Every kept entry gets a row even with no Pfam xref, because the map's
        # name -> accession half is what resolves HIPPIE identifiers.
        assert pfam[2] == "Q11111\tNOGO_HUMAN\tb0001 CG42\t"

        fasta = got["seq.fasta.gz"].splitlines()
        assert fasta[0].startswith(">sp|O00192|ARVC_HUMAN ")
        assert fasta[0].endswith(" OS=Homo sapiens (Human) OX=9606")
        assert fasta[1] == "MSSRSVSRSRARSRSPRRSRSRSRSPAYKR"

        # The mouse entry is gone from all three.
        for name, text in got.items():
            assert "P99999" not in text and "MOUSY" not in text, name


def test_no_filter_keeps_every_species():
    with tempfile.TemporaryDirectory() as tmp:
        got = run(tmp, "")
        assert "P99999\tGO:0000001" in got["go.tsv.gz"]
        assert "P99999\tMOUSY_MOUSE\tMus1\tPF99999;" in got["pfam.tsv.gz"]
        assert ">sp|P99999|MOUSY_MOUSE" in got["seq.fasta.gz"]


def test_taxon_matching_nothing_is_fatal():
    """A silently empty parse would hand the enrichment chain three empty files
    and produce a DB with no GO terms and no sequences, rather than an error."""
    with tempfile.TemporaryDirectory() as tmp:
        dat = os.path.join(tmp, "sprot.dat.gz")
        with gzip.open(dat, "wt") as fh:
            fh.write(DAT)
        proc = subprocess.run(
            [sys.executable, PARSER, "--dat", dat,
             "--go-terms", os.path.join(tmp, "go.tsv.gz"),
             "--pfam-map", os.path.join(tmp, "pfam.tsv.gz"),
             "--sequences", os.path.join(tmp, "seq.fasta.gz"),
             "--string-map", os.path.join(tmp, "string.tsv.gz"),
             "--taxon-ids", "7227",
             "--versions", os.path.join(tmp, "versions.yml"),
             "--process-name", "TEST"],
            capture_output=True, text=True,
        )
        assert proc.returncode != 0
        assert "no entries survived" in proc.stderr, proc.stderr


if __name__ == "__main__":
    test_human_filter_and_output_shapes()
    test_string_map_replaces_the_idmapping_download()
    test_no_filter_keeps_every_species()
    test_taxon_matching_nothing_is_fatal()
    print("OK: parse_swissprot_dat emits the four shapes its consumers parse")
