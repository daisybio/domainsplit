#!/usr/bin/env python3
"""Parse the UniProt flat files that form the run's protein universe into the
four TSV/FASTA inputs the enrichment chain used to fetch separately.

One pass over static FTP files (``uniprot_sprot_human.dat.gz``,
``uniprot_sprot.dat.gz``, ``uniprot_trembl_human.dat.gz``, ... -- whichever
``--instance_tier`` resolved to) replaces two ``rest.uniprot.org`` stream
queries, which are throttled to ~4 KB/s per connection and drop the connection
long before an all-reviewed query finishes. The outputs are byte-compatible with
what the old sources produced, so ``insert_protein_go_terms.py``,
``build_swissprot_pfam_map.py`` and ``insert_protein_sequences.py`` are unchanged:

* ``--go-terms``   ``Entry\\tGene Ontology IDs`` TSV.gz, GO ids joined by "; ".
                   Every ``DR   GO;`` line is kept, evidence code included, which
                   is what the REST ``go_id`` field returned. **All** entries:
                   TrEMBL GO is IEA, the DB records it as protein annotation and
                   nothing splits on it.
* ``--pfam-map``   ``Entry\\tEntry Name\\tGene Names\\tPfam`` TSV.gz, Pfam
                   accessions joined by ";" with a trailing ";", gene names
                   space-separated -- the shape of the REST
                   ``accession,id,gene_names,xref_pfam`` stream. **Reviewed
                   entries only** by default; see ``--pfam-map-reviewed-only``.
* ``--sequences``  ``>sp|ACC|ENTRY_NAME ...`` / ``>tr|ACC|ENTRY_NAME ...``
                   FASTA.gz, the same headers UniProt's own FASTA carries. The
                   prefix follows the entry's review flag: ``insert_protein_
                   sequences.py`` splits on "|" and tolerates either, but a
                   TrEMBL entry labelled ``sp|`` is a lie in a published file.
* ``--string-map`` ``Entry\tSTRING`` TSV.gz, one row per ``DR   STRING;`` id, in
                   the ``uniprot_id\tid_type\tsymbol`` shape ``insert_ppi.py``
                   already reads. STRING ids carry their own taxon prefix
                   (``9606.ENSP…``) and are globally unique, so a multi-species
                   map joins correctly with no taxon bookkeeping -- and
                   FETCH_STRING_LINKS reads the taxon list straight off those
                   prefixes.

``--dat`` is repeatable and the files are parsed in the given order. **Swiss-Prot
must be listed first**: an accession promoted from TrEMBL to Swiss-Prot between
releases exists in both files, and the reviewed record is the one to keep, so a
duplicate accession is resolved first-writer-wins rather than by whichever file
happened to be read last.

``--taxon-ids`` restricts the output to the given NCBI taxonomy ids; omit it (or
pass an empty string) to keep every entry. Both it and the file list are derived
from ``--instance_tier``, which is the single knob that decides the protein
universe for the whole run: the same list goes to FETCH_DOMAIN_META, and
``ingest_instances.py`` asserts against the same taxa. Nothing here is
human-specific.
"""

import argparse
import gzip
import re
import sys
from collections import Counter

# "OX   NCBI_TaxID=9606;" -- an evidence tag may follow the id.
TAXID_RE = re.compile(r"NCBI_TaxID=(\d+)")
# UniProt evidence tags: "Name=TP53 {ECO:0000255|HAMAP-Rule:MF_00001};"
EVIDENCE_RE = re.compile(r"\s*\{[^}]*\}")

GENE_KEYS = ("Name=", "Synonyms=", "OrderedLocusNames=", "ORFNames=")

#: The ID line's second whitespace field. "ID   TP53_HUMAN   Reviewed;   393 AA."
REVIEW_TOKENS = {"Reviewed;": True, "Unreviewed;": False}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dat", required=True, action="append", default=[],
                   help="UniProt flat file (.dat.gz); repeatable, Swiss-Prot first")
    p.add_argument("--go-terms", required=True, help="output TSV.gz")
    p.add_argument("--pfam-map", required=True, help="output TSV.gz")
    p.add_argument("--sequences", required=True, help="output FASTA.gz")
    p.add_argument("--string-map", required=True, help="output TSV.gz: Entry -> STRING id")
    p.add_argument("--taxon-ids", default="",
                   help="comma-separated NCBI taxon ids; empty keeps every entry")
    p.add_argument(
        "--pfam-map-reviewed-only",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="write --pfam-map for reviewed entries only (default). That output feeds "
        "build_swissprot_pfam_map.py -> parse_single_domain_ppi.py, which infers a DDI "
        "from a HIPPIE PPI between two *single-domain* proteins. On a TrEMBL entry the "
        "Pfam annotation is automatic and the sequence is often a fragment, so "
        "'this protein has exactly one domain' is not a claim that source can support "
        "-- it would manufacture DDIs. Unreviewed proteins are instance "
        "representatives, not evidence.",
    )
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


def flush(entry, out, stats, seen, pfam_map_reviewed_only):
    """Write one finished entry to the four outputs.

    ``seen`` carries every accession written so far, across all input files:
    first writer wins, which is why Swiss-Prot has to be listed first.
    """
    acc = entry["acc"]
    if not acc:
        stats["no_accession"] += 1
        return
    if acc in seen:
        stats["duplicate_accession"] += 1
        return
    seen.add(acc)

    reviewed = entry["reviewed"]
    flag = "reviewed" if reviewed else "unreviewed"
    stats["kept"] += 1
    stats[f"kept_{flag}"] += 1

    if entry["go"]:
        out["go"].write(f"{acc}\t{'; '.join(entry['go'])}\n")
        stats["with_go"] += 1

    if reviewed or not pfam_map_reviewed_only:
        pfams = ";".join(entry["pfam"])
        out["pfam"].write(
            f"{acc}\t{entry['id']}\t{' '.join(gene_tokens(entry['gn']))}\t"
            f"{pfams + ';' if pfams else ''}\n"
        )
        stats[f"pfam_map_{flag}"] += 1
        if len(entry["pfam"]) == 1:
            stats["single_pfam"] += 1
    else:
        stats["pfam_map_skipped_unreviewed"] += 1

    seq = "".join(entry["seq"])
    if seq:
        organism = entry["os"].rstrip(".")
        prefix = "sp" if reviewed else "tr"
        out["fasta"].write(
            f">{prefix}|{acc}|{entry['id']} {entry['de']} OS={organism} OX={entry['ox']}\n"
        )
        for i in range(0, len(seq), 60):
            out["fasta"].write(seq[i:i + 60] + "\n")
        stats["with_sequence"] += 1

    # One row per STRING id, matching load_uniprot_id_mapping()'s
    # "uniprot \t id_type \t symbol" expectation. An entry usually has one; a few
    # carry several and every one of them is a valid join key. Unreviewed entries
    # almost never have any -- UniProt does not cross-reference STRING for them --
    # which is why the reviewed/unreviewed split is reported below.
    for string_id in entry["string"]:
        out["string"].write(f"{acc}\tSTRING\t{string_id}\n")
    if entry["string"]:
        stats["with_string"] += 1
        stats[f"with_string_{flag}"] += 1


def new_entry():
    return {"acc": "", "id": "", "de": "", "os": "", "ox": "", "reviewed": None,
            "gn": [], "go": [], "pfam": [], "string": [], "seq": [], "in_seq": False}


def review_flag(line, path):
    """``True``/``False`` off an ``ID`` line; fatal when the token is missing.

    Not defaulted to ``False``: a silently mis-parsed flag would put every TrEMBL
    protein in the reviewed stratum, and there would be no symptom until someone
    looked at ``protein.reviewed``.
    """
    fields = line[5:].split()
    if len(fields) >= 2 and fields[1] in REVIEW_TOKENS:
        return REVIEW_TOKENS[fields[1]]
    raise SystemExit(
        f"ERROR: {path}: ID line carries neither 'Reviewed;' nor 'Unreviewed;', so the "
        f"entry's review status cannot be read: {line.rstrip()!r}"
    )


def parse_file(path, wanted, out, stats, seen, pfam_map_reviewed_only):
    """One streaming pass over one flat file. Returns its own entry counter."""
    per_file = Counter()
    entry = new_entry()
    # Set once OX is read and the entry is unwanted: every remaining line of that
    # entry is then skipped with a single two-byte compare, which is what makes a
    # filtered pass over hundreds of thousands of entries cheap.
    skipping = False

    with gzip.open(path, "rt", encoding="latin-1") as dat:
        for line in dat:
            tag = line[:2]

            if tag == "//":
                stats["entries"] += 1
                per_file["entries"] += 1
                if not skipping:
                    before = stats["kept"]
                    flush(entry, out, stats, seen, pfam_map_reviewed_only)
                    per_file["kept"] += stats["kept"] - before
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
                entry["reviewed"] = review_flag(line, path)
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
            stats["truncated_entry_dropped"] += 1

    return per_file


def main():
    args = parse_args()
    wanted = {t.strip() for t in args.taxon_ids.split(",") if t.strip()}
    print(
        f"[parse_swissprot] {len(args.dat)} flat file(s): {', '.join(args.dat)}",
        flush=True,
    )
    print(
        f"[parse_swissprot] taxon filter: {sorted(wanted) if wanted else 'none (all species)'}"
        f"; pfam map: {'reviewed entries only' if args.pfam_map_reviewed_only else 'every entry'}",
        flush=True,
    )

    stats = Counter()
    seen = set()

    with (
        gzip.open(args.go_terms, "wt", newline="") as out_go,
        gzip.open(args.pfam_map, "wt", newline="") as out_pfam,
        gzip.open(args.sequences, "wt", newline="") as out_fasta,
        gzip.open(args.string_map, "wt", newline="") as out_string,
    ):
        out = {"go": out_go, "pfam": out_pfam, "fasta": out_fasta, "string": out_string}
        out_go.write("Entry\tGene Ontology IDs\n")
        out_pfam.write("Entry\tEntry Name\tGene Names\tPfam\n")

        for path in args.dat:
            per_file = parse_file(path, wanted, out, stats, seen, args.pfam_map_reviewed_only)
            print(
                f"[parse_swissprot] {path}: {per_file['entries']} entries, "
                f"{per_file['kept']} kept",
                flush=True,
            )
            # A file that contributes nothing is a configuration error worth
            # naming: the wrong division file, or a taxon filter that excludes
            # everything it holds.
            if per_file["kept"] == 0:
                raise SystemExit(
                    f"ERROR: {path} contributed no entry. Either it is not the file "
                    f"--instance_tier meant to use, or --taxon-ids '{args.taxon_ids}' "
                    "excludes everything in it (a second flat file whose accessions are "
                    "all already in the first would also do this)."
                )

    print(
        "[parse_swissprot] " + " ".join(f"{k}={v}" for k, v in sorted(stats.items())),
        flush=True,
    )
    if stats["kept"] == 0:
        raise SystemExit(
            "ERROR: no entries survived parsing. Either the .dat files are truncated "
            f"or --taxon-ids '{args.taxon_ids}' matches nothing in them."
        )
    if stats["kept_unreviewed"] and not stats["with_string_unreviewed"]:
        print(
            f"[parse_swissprot] {stats['kept_unreviewed']} unreviewed entries carry no "
            "`DR STRING;` line at all -- UniProt does not cross-reference STRING for "
            "them, so those proteins get no PPI enrichment and no links file changes "
            "that",
            flush=True,
        )

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
