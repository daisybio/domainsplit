#!/usr/bin/env python3
"""Parse the UniProt-SwissProt flat file into the four TSV/FASTA inputs the
enrichment chain used to fetch separately.

One pass over ``uniprot_sprot.dat.gz`` (a static FTP file that transfers at full
line speed) replaces two ``rest.uniprot.org`` stream queries, which are throttled
to ~4 KB/s per connection and drop the connection long before an all-reviewed
query finishes. The three outputs are byte-compatible with what the old sources
produced, so ``insert_protein_go_terms.py``, ``build_swissprot_pfam_map.py`` and
``insert_protein_sequences.py`` are unchanged:

* ``--go-terms``   ``Entry\\tGene Ontology IDs`` TSV.gz, GO ids joined by "; ".
                   Every ``DR   GO;`` line is kept, evidence code included, which
                   is what the REST ``go_id`` field returned.
* ``--pfam-map``   ``Entry\\tEntry Name\\tGene Names\\tPfam`` TSV.gz, Pfam
                   accessions joined by ";" with a trailing ";", gene names
                   space-separated -- the shape of the REST
                   ``accession,id,gene_names,xref_pfam`` stream.
* ``--sequences``  ``>sp|ACC|ENTRY_NAME ...`` FASTA.gz, the same headers
                   ``uniprot_sprot.fasta.gz`` carries.
* ``--string-map`` ``Entry\tSTRING`` TSV.gz, one row per ``DR   STRING;`` id, in
                   the ``uniprot_id\tid_type\tsymbol`` shape ``insert_ppi.py``
                   already reads -- so it is a drop-in replacement for the
                   per-organism ``<ORG>_<taxid>_idmapping.dat.gz`` download, and
                   unlike that file it covers **every** reviewed species from the
                   one pass we were already making. STRING ids carry their own
                   taxon prefix (``9606.ENSP…``) and are globally unique, so a
                   multi-species map joins correctly with no taxon bookkeeping.

``--taxon-ids`` restricts the output to the given NCBI taxonomy ids; omit it (or
pass an empty string) to keep every reviewed entry. The default at the pipeline
level is 9606, and it is the single knob that decides the protein universe for the
whole run: ``instance_tier`` selects from within it, and ``ingest_instances.py``
asserts against it. Nothing here is human-specific.
"""

import argparse
import gzip
import re
import sys

# "OX   NCBI_TaxID=9606;" -- an evidence tag may follow the id.
TAXID_RE = re.compile(r"NCBI_TaxID=(\d+)")
# UniProt evidence tags: "Name=TP53 {ECO:0000255|HAMAP-Rule:MF_00001};"
EVIDENCE_RE = re.compile(r"\s*\{[^}]*\}")

GENE_KEYS = ("Name=", "Synonyms=", "OrderedLocusNames=", "ORFNames=")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dat", required=True, help="uniprot_sprot.dat.gz")
    p.add_argument("--go-terms", required=True, help="output TSV.gz")
    p.add_argument("--pfam-map", required=True, help="output TSV.gz")
    p.add_argument("--sequences", required=True, help="output FASTA.gz")
    p.add_argument("--string-map", required=True, help="output TSV.gz: Entry -> STRING id")
    p.add_argument("--taxon-ids", default="",
                   help="comma-separated NCBI taxon ids; empty keeps every entry")
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def gene_tokens(gn_lines):
    """The REST ``gene_names`` field: name, synonyms, ordered locus and ORF names
    of every gene on the entry, whitespace-separated, evidence tags stripped."""
    tokens = []
    for raw in gn_lines:
        for chunk in EVIDENCE_RE.sub("", raw).split(";"):
            chunk = chunk.strip()
            if not chunk:
                continue
            for key in GENE_KEYS:
                if chunk.startswith(key):
                    for token in chunk[len(key):].split(","):
                        token = token.strip()
                        if token and token not in tokens:
                            tokens.append(token)
                    break
    return tokens


def flush(entry, out_go, out_pfam, out_fasta, out_string, stats):
    """Write one finished entry to the four outputs."""
    acc = entry["acc"]
    if not acc:
        stats["no_accession"] += 1
        return
    stats["kept"] += 1

    if entry["go"]:
        out_go.write(f"{acc}\t{'; '.join(entry['go'])}\n")
        stats["with_go"] += 1

    pfams = ";".join(entry["pfam"])
    out_pfam.write(
        f"{acc}\t{entry['id']}\t{' '.join(gene_tokens(entry['gn']))}\t"
        f"{pfams + ';' if pfams else ''}\n"
    )
    if len(entry["pfam"]) == 1:
        stats["single_pfam"] += 1

    seq = "".join(entry["seq"])
    if seq:
        organism = entry["os"].rstrip(".")
        out_fasta.write(f">sp|{acc}|{entry['id']} {entry['de']} OS={organism} OX={entry['ox']}\n")
        for i in range(0, len(seq), 60):
            out_fasta.write(seq[i:i + 60] + "\n")
        stats["with_sequence"] += 1

    # One row per STRING id, matching load_uniprot_id_mapping()'s
    # "uniprot \t id_type \t symbol" expectation. An entry usually has one; a few
    # carry several and every one of them is a valid join key.
    for string_id in entry["string"]:
        out_string.write(f"{acc}\tSTRING\t{string_id}\n")
    if entry["string"]:
        stats["with_string"] += 1


def new_entry():
    return {"acc": "", "id": "", "de": "", "os": "", "ox": "",
            "gn": [], "go": [], "pfam": [], "string": [], "seq": [], "in_seq": False}


def main():
    args = parse_args()
    wanted = {t.strip() for t in args.taxon_ids.split(",") if t.strip()}
    print(
        f"[parse_swissprot] taxon filter: {sorted(wanted) if wanted else 'none (all species)'}",
        flush=True,
    )

    stats = dict.fromkeys(
        ("entries", "kept", "with_go", "with_sequence", "with_string", "single_pfam",
         "wrong_taxon", "no_accession"), 0
    )
    entry = new_entry()
    # Set once OX is read and the entry is unwanted: every remaining line of that
    # entry is then skipped with a single two-byte compare, which is what makes a
    # filtered pass over ~570k entries cheap.
    skipping = False

    with (
        gzip.open(args.dat, "rt", encoding="latin-1") as dat,
        gzip.open(args.go_terms, "wt", newline="") as out_go,
        gzip.open(args.pfam_map, "wt", newline="") as out_pfam,
        gzip.open(args.sequences, "wt", newline="") as out_fasta,
        gzip.open(args.string_map, "wt", newline="") as out_string,
    ):
        out_go.write("Entry\tGene Ontology IDs\n")
        out_pfam.write("Entry\tEntry Name\tGene Names\tPfam\n")

        for line in dat:
            tag = line[:2]

            if tag == "//":
                stats["entries"] += 1
                if not skipping:
                    flush(entry, out_go, out_pfam, out_fasta, out_string, stats)
                entry = new_entry()
                skipping = False
                continue

            if skipping:
                continue

            if entry["in_seq"]:
                entry["seq"].append(line.strip().replace(" ", ""))
                continue

            if tag == "ID":
                entry["id"] = line[5:].split()[0]
            elif tag == "AC":
                if not entry["acc"]:
                    # The first accession on the first AC line is the primary one;
                    # the rest are secondary and are what the REST TSV omits.
                    entry["acc"] = line[5:].split(";")[0].strip()
            elif tag == "DE":
                if not entry["de"]:
                    body = line[5:].strip()
                    if body.startswith(("RecName: Full=", "SubName: Full=")):
                        entry["de"] = EVIDENCE_RE.sub("", body.split("Full=", 1)[1]).rstrip(";")
            elif tag == "GN":
                entry["gn"].append(line[5:].strip())
            elif tag == "OS":
                entry["os"] += (" " if entry["os"] else "") + line[5:].strip()
            elif tag == "OX":
                m = TAXID_RE.search(line)
                entry["ox"] = m.group(1) if m else ""
                if wanted and entry["ox"] not in wanted:
                    stats["wrong_taxon"] += 1
                    skipping = True
            elif tag == "DR":
                body = line[5:]
                if body.startswith("GO; "):
                    entry["go"].append(body[4:].split(";")[0].strip())
                elif body.startswith("Pfam; "):
                    pfam = body[6:].split(";")[0].strip()
                    if pfam and pfam not in entry["pfam"]:
                        entry["pfam"].append(pfam)
                elif body.startswith("STRING; "):
                    # "DR   STRING; 9606.ENSP00000269305; -."
                    string_id = body[8:].split(";")[0].strip()
                    if string_id and string_id not in entry["string"]:
                        entry["string"].append(string_id)
            elif tag == "SQ":
                entry["in_seq"] = True

        # A file truncated mid-entry has no closing "//"; that entry is dropped
        # rather than written half-parsed.
        if entry["acc"] and not skipping:
            stats["truncated_entry_dropped"] = 1

    print(
        "[parse_swissprot] "
        + " ".join(f"{k}={v}" for k, v in stats.items()),
        flush=True,
    )
    if stats["kept"] == 0:
        raise SystemExit(
            "ERROR: no entries survived parsing. Either the .dat file is truncated "
            f"or --taxon-ids '{args.taxon_ids}' matches nothing in it."
        )

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
