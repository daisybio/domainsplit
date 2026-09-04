#!/usr/bin/env python3
"""Cut the `-profile test` fixtures out of the real source files.

The real sources are multi-GB and not in git; this script reduces them to a
handful of Pfam families and writes everything `-profile test` needs into
``tests/data/``. Run it again when a source is refreshed:

    tests/bin/make_test_fixtures.py --sources ~/Downloads/pipeline_files

Provenance of the inputs is documented in ``tests/data/README.md``.

## How the families are chosen

Everything keys off one decision: which Pfam families the fixture covers. A
family is only useful if it has **human** domain instances, because the pipeline
runs `instance_tier = human_reviewed` and a family with no human instance is dropped
by PRUNE_UNREPRESENTED_DDIS. So the selection is:

1. build the 3did DDI graph over Pfam accessions (only 3did's `DDI1` and
   `Domain` tables are read by the pipeline, so only those are subset);
2. shortlist families by 3did degree, preferring ones PPIDM or Negatome also
   mention -- that is what gives the external test set something to hold;
3. take **one** pass over Pfam-A.fasta.gz counting human records per shortlisted
   family, and keep the families that clear `--min-human-instances`;
4. induce the 3did subgraph on the survivors.

Step 3 is the expensive part (~6.3 GB of gzip) and the reason this is a script
rather than a handful of `grep` invocations.

## What is real and what is synthesised

Real subsets, cut from the files under ``--sources``: the 3did dump, HIPPIE,
PPIDM, Negatome, the Y2H/MS parquet, Pfam-A.fasta.gz, Pfam-A.clans.tsv.

Synthesised, because no local copy of the source exists and the pipeline only
needs the *format*: the reviewed-human UniProt TSV, `speclist.txt`,
`reviewed.list`, and the gene → UniProt → Pfam mapping that lets
BUILD_CANDIDATE_NETWORK skip the UniProt API. They are all derived from the real
records above, so they stay mutually consistent: the UniProt TSV's accessions are
the accessions Pfam-A.fasta actually carries, and the mapping's gene names are the
gene names the parquet actually contains.
"""

import argparse
import csv
import gzip
import io
import json
import os
import random
import re
import sys
from collections import Counter, defaultdict

import universe_tier_fixtures

# Pfam-A.fasta header: "A0A0A0MRZ7_HUMAN/1-30 A0A0A0MRZ7.1 PF10417.14;..."
# Same shape ppi-splitting's fetch_domains.py parses, so records copied verbatim
# out of the real file are guaranteed to survive its reader.
PFAM_HEADER_RE = re.compile(
    r"^(?P<entry>\S+)/(?P<start>\d+)-(?P<end>\d+)\s+(?P<acc>[A-Za-z0-9]+)(?:\.\d+)?\s+(?P<family>PF\d+)"
)

# 3did tables insert_3did.py reads. Everything else in the 1.3 GB dump is dropped.
THREEDID_TABLES = ("DDI1", "Domain")


def log(msg):
    print(f"[fixtures] {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sources", required=True, help="directory holding the real source files")
    p.add_argument("--out", default=None, help="fixture directory (default: <repo>/tests/data)")
    p.add_argument("--families", type=int, default=85,
                   help="number of Pfam families to keep. Raised from 40 so `ilp_candidates` is "
                        "usable at fixture scale: candidate pairs survive a split only if *both* "
                        "families land in it, so survival is quadratic in the split's share and a "
                        "0.1 val split of 40 families could not cover its own positives. 85 rather "
                        "than the ~60 that first clears the bar, because the margin is what matters: "
                        "measured over the fixture's own pair set, P(candidates < positives) for a "
                        "0.1 split is 0.49% at 70 families and 0.02% at 85, for 21% more DDIs -- and "
                        "the real ILP partition is lumpier than the uniform draw that measures")
    p.add_argument("--shortlist", type=int, default=300,
                   help="families carried into the Pfam-A.fasta pass, before the human filter")
    p.add_argument("--min-human-instances", type=int, default=3,
                   help="human records a family needs to be kept; 2 is the minimum for a "
                        "self-pair DDI to have an instance pair with distinct parents")
    p.add_argument("--max-records-per-family", type=int, default=40,
                   help="records copied per family into the miniature Pfam-A.fasta.gz")
    p.add_argument("--hippie-pairs", type=int, default=40)
    p.add_argument("--parquet-rows", type=int, default=400)
    p.add_argument("--pfam-release", default="38.0",
                   help="release string; only labels the interpro_cache directory")
    p.add_argument("--seed", type=int, default=42)
    return p.parse_args()


# ---------------------------------------------------------------------------
# 3did
# ---------------------------------------------------------------------------


def split_sql_tuples(values):
    """Yield the field lists of ``('a','b',NULL),('c',...)`` without eval()."""
    fields, current, depth = [], [], 0
    buf = io.StringIO()
    in_str = False
    i = 0
    while i < len(values):
        ch = values[i]
        if in_str:
            if ch == "\\":
                buf.write(values[i : i + 2])
                i += 2
                continue
            if ch == "'":
                in_str = False
            else:
                buf.write(ch)
            i += 1
            continue
        if ch == "'":
            in_str = True
        elif ch == "(":
            depth += 1
            if depth == 1:
                current = []
                buf = io.StringIO()
        elif ch == ")":
            current.append(buf.getvalue())
            buf = io.StringIO()
            depth = 0
            fields.append(current)
        elif ch == "," and depth == 1:
            current.append(buf.getvalue())
            buf = io.StringIO()
        elif depth == 1:
            buf.write(ch)
        i += 1
    return fields


def read_3did(path):
    """``(preamble, {table: create_block}, ddi_name_pairs, domain_rows)``.

    The CREATE blocks are copied verbatim so the fixture stays a real mysqldump
    that ``bin/mysql2sqlite`` handles exactly like the full one.
    """
    preamble, blocks = [], {}
    ddi_pairs, domain_rows = [], []
    table = None
    with open(path, encoding="utf8", errors="replace") as fh:
        for line in fh:
            m = re.match(r"^CREATE TABLE `(\w+)`", line)
            if m:
                table = m.group(1)
                if table in THREEDID_TABLES:
                    blocks[table] = [line]
                continue
            if table is None:
                preamble.append(line)
                continue
            if table in THREEDID_TABLES and not line.startswith("INSERT INTO"):
                # keep the structural tail (column list, engine line, SET pragmas)
                # but stop before the next table's marker
                if not re.match(r"^(--|LOCK TABLES|UNLOCK TABLES|/\*!40000)", line):
                    blocks[table].append(line)
                continue
            if line.startswith("INSERT INTO"):
                m = re.match(r"^INSERT INTO `(\w+)` VALUES (.*);\s*$", line)
                if not m or m.group(1) not in THREEDID_TABLES:
                    continue
                rows = split_sql_tuples(m.group(2))
                if m.group(1) == "DDI1":
                    ddi_pairs.extend((r[0], r[2]) for r in rows if len(r) > 2)
                else:
                    domain_rows.extend(rows)
    log(f"3did: {len(ddi_pairs)} DDI rows, {len(domain_rows)} domain rows")
    return preamble, blocks, ddi_pairs, domain_rows


def write_3did_fixture(out, preamble, blocks, keep_ddi, keep_domain):
    def insert(table, rows):
        body = ",".join(
            "(" + ",".join("NULL" if f == "NULL" else "'" + f.replace("'", "''") + "'" for f in row) + ")"
            for row in rows
        )
        return f"INSERT INTO `{table}` VALUES {body};\n"

    path = os.path.join(out, "3did.sql.gz")
    with gzip.open(path, "wt", encoding="utf8") as fh:
        fh.writelines(preamble)
        for table, rows in (("DDI1", keep_ddi), ("Domain", keep_domain)):
            fh.write(f"\n--\n-- Table structure for table `{table}`\n--\n\n")
            fh.write(f"DROP TABLE IF EXISTS `{table}`;\n")
            fh.writelines(blocks[table])
            fh.write(f"\n--\n-- Dumping data for table `{table}`\n--\n\n")
            fh.write(f"LOCK TABLES `{table}` WRITE;\n")
            if rows:
                fh.write(insert(table, rows))
            fh.write("UNLOCK TABLES;\n")
    log(f"wrote {path} ({os.path.getsize(path)} bytes)")


# ---------------------------------------------------------------------------
# external sources
# ---------------------------------------------------------------------------


def pfam_of(token):
    """``PF00001`` out of ``PF00001.12``, ``10114/PF00069`` or ``PF00001``, else None.

    Same normalisation the parsers use: PPIDM prefixes each accession with an
    Entrez id, 3did suffixes a version.
    """
    stem = token.strip().split("/")[-1].split(".")[0]
    return stem if re.fullmatch(r"PF\d{5}", stem) else None


def read_ppidm(path):
    """``(header, {(pfam_a, pfam_b): {class, ...}})`` from predicted_ddi_ppi.tsv.

    A pair predicted under two classes keeps both, so the fixture can cover
    INSERT_EXTERNAL_SOURCES' source merging (`PPIDM,PPIDM_Gold,PPIDM_Silver`).
    """
    pairs = defaultdict(set)
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) < 3:
                continue
            a, b = pfam_of(cols[0]), pfam_of(cols[1])
            if a and b:
                pairs[(a, b)].add(cols[2].strip())
    log(f"ppidm: {len(pairs)} pairs, "
        f"{sum(1 for c in pairs.values() if len(c) > 1)} with more than one class")
    return header, pairs


def read_negatome(path):
    pairs = []
    with open(path) as fh:
        for line in fh:
            tokens = line.split()
            if len(tokens) < 2:
                continue
            a, b = pfam_of(tokens[0]), pfam_of(tokens[1])
            if a and b:
                pairs.append((a, b))
    log(f"negatome: {len(pairs)} pairs")
    return pairs


# ---------------------------------------------------------------------------
# Pfam-A.fasta
# ---------------------------------------------------------------------------


def scan_pfam_fasta(path, wanted, max_per_family):
    """One pass, collecting up to ``max_per_family`` human records per family.

    Returns ``({family: [(header, sequence)]}, {family: n_human_seen})``. Human
    records only: the pipeline runs `instance_tier = human_reviewed`, so a non-human
    record in the fixture would be dead weight in a 6.3 GB pass that is already
    the slow step.
    """
    kept = defaultdict(list)
    seen = Counter()
    with gzip.open(path, "rt", encoding="ascii", errors="replace") as fh:
        header, chunks, family = None, [], None
        for line in fh:
            if line.startswith(">"):
                if header and family and len(kept[family]) < max_per_family:
                    kept[family].append((header, "".join(chunks)))
                header, chunks, family = None, [], None
                m = PFAM_HEADER_RE.match(line[1:].rstrip())
                if not m or m["family"] not in wanted:
                    continue
                if not m["entry"].endswith("_HUMAN"):
                    continue
                seen[m["family"]] += 1
                header, family = line.rstrip(), m["family"]
            elif header:
                chunks.append(line.strip())
        if header and family and len(kept[family]) < max_per_family:
            kept[family].append((header, "".join(chunks)))
    log(f"pfam fasta: {sum(seen.values())} human records over {len(seen)} shortlisted families")
    return kept, seen


def write_pfam_fixtures(out, records, clans_path, families, release):
    pfam_dir = os.path.join(out, "pfam")
    os.makedirs(pfam_dir, exist_ok=True)

    fasta = os.path.join(pfam_dir, "Pfam-A.fasta.gz")
    with gzip.open(fasta, "wt", encoding="ascii") as fh:
        # Family blocks contiguous, as in the real file.
        for family in families:
            for header, seq in records[family]:
                fh.write(f"{header}\n")
                for i in range(0, len(seq), 60):
                    fh.write(seq[i : i + 60] + "\n")
    log(f"wrote {fasta} ({os.path.getsize(fasta)} bytes)")

    clans = os.path.join(pfam_dir, "Pfam-A.clans.tsv.gz")
    wanted = set(families)
    n = 0
    with open(clans_path) as src, gzip.open(clans, "wt", encoding="utf8") as fh:
        for line in src:
            if pfam_of(line.split("\t")[0] or "") in wanted:
                fh.write(line)
                n += 1
    log(f"wrote {clans} ({n} of {len(wanted)} families have a clans row)")
    return pfam_dir


def write_interpro_cache(pfam_dir, release, mnemonics, accessions):
    """Pre-populate ppi-splitting's download cache so FETCH_DOMAIN_META runs offline.

    `fetch_domains.py` reads `speclist.txt` and `reviewed.list` through a
    release-keyed cache and only downloads a file the cache does not have, so
    writing them here removes the two remaining network calls (the third,
    Pfam.version, is skipped by pinning `--pfam_release`).
    """
    cache = os.path.join(pfam_dir, "interpro_cache", f"pfam-{release}")
    os.makedirs(cache, exist_ok=True)
    with open(os.path.join(cache, "speclist.txt"), "w") as fh:
        fh.write("Real code Taxon Node/Scientific name\n")
        for mnemonic in sorted(mnemonics):
            taxon = 9606 if mnemonic == "HUMAN" else 32644
            fh.write(f"{mnemonic} E {taxon}: N=Homo sapiens\n")
    with open(os.path.join(cache, "reviewed.list"), "w") as fh:
        # Every accession in the fixture counts as reviewed, so each instance
        # lands in the `human_reviewed` stratum -- the one every tier fills first.
        for acc in sorted(accessions):
            fh.write(f"{acc}\n")
    # Pfam-A.dead is fetched when an accession fails to resolve. Empty is a valid
    # dead list (the parser looks for `#=GF AC` lines) and closes the last hole
    # through which a test run could reach the network.
    with open(os.path.join(cache, "Pfam-A.dead"), "w") as fh:
        fh.write("")
    log(f"wrote {cache}/{{speclist.txt,reviewed.list,Pfam-A.dead}} "
        f"({len(accessions)} reviewed accessions)")
    return cache


# ---------------------------------------------------------------------------
# HIPPIE + the reviewed-human UniProt TSV
# ---------------------------------------------------------------------------


def build_hippie_and_swissprot(hippie_path, families, n_pairs, min_score, rng):
    """Real HIPPIE rows, plus the synthetic UniProt TSV that makes them resolve.

    PARSE_SINGLE_DOMAIN_PPI keeps a HIPPIE row only when *both* identifiers are
    single-domain proteins, so the two fixtures have to agree. Rather than hunt
    for real rows whose real proteins happen to carry our families, the HIPPIE
    rows stay real and the UniProt TSV is written to describe exactly their
    identifiers, each as a single-domain protein of one of our families.
    """
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))), "bin"))
    from parse_single_domain_ppi import detect_layout

    rows, ids = [], []
    with open(hippie_path) as fh:
        first = fh.readline()
        (col_a, col_b, col_score), is_header = detect_layout(first)
        header = first if is_header else None
        lines = fh if is_header else [first] + list(fh)
        for line in lines:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= col_score:
                continue
            try:
                if float(cols[col_score]) < min_score:
                    continue
            except ValueError:
                continue
            rows.append(line)
            # One protein per side, as (accession, entry name) -- the two columns
            # of the current layout describe the same protein, so they must not
            # become two synthetic entries with two different families.
            for side in (col_a, col_b):
                tokens = [cols[i] for i in side if cols[i]]
                if tokens:
                    ids.append((tokens[0], tokens[-1]))
            if len(rows) >= n_pairs:
                break
    if header:
        rows.insert(0, header)
    unique_ids = sorted(set(ids))
    log(f"hippie: {len(rows)} rows over {len(unique_ids)} proteins")

    # Exactly one family per protein, so every kept HIPPIE row yields a DDI.
    # Families cycle, so the rows produce several distinct pairs rather than one.
    tsv_rows = []
    for i, (acc, name) in enumerate(unique_ids):
        family = families[i % len(families)]
        tsv_rows.append((acc, name, name.split("_")[0], family))
    return rows, tsv_rows


def write_hippie(out, rows):
    path = os.path.join(out, "hippie.tsv")
    with open(path, "w") as fh:
        fh.writelines(rows)
    log(f"wrote {path} ({len(rows)} rows)")


def write_dat_entry(fh, acc, entry_name, gene, seq, pfams=(), go=(), taxid="9606"):
    """One UniProt-SwissProt flat-file entry, in the subset of the format
    `parse_swissprot_dat.py` reads: ID/AC/DE/GN/OS/OX, the DR GO and DR Pfam
    cross-references, and the SQ block."""
    fh.write(f"ID   {entry_name}   Reviewed;   {len(seq)} AA.\n")
    fh.write(f"AC   {acc};\n")
    fh.write(f"DE   RecName: Full=Synthetic fixture protein {acc};\n")
    if gene:
        fh.write(f"GN   Name={gene};\n")
    fh.write("OS   Homo sapiens (Human).\n")
    fh.write(f"OX   NCBI_TaxID={taxid};\n")
    for go_id in go:
        fh.write(f"DR   GO; {go_id}; F:synthetic fixture term; IEA:fixture.\n")
    for pfam in pfams:
        fh.write(f"DR   Pfam; {pfam}; {pfam}_name; 1.\n")
    fh.write(f"SQ   SEQUENCE   {len(seq)} AA;  0 MW;  0000000000000000 CRC64;\n")
    for i in range(0, len(seq), 60):
        fh.write("     " + seq[i : i + 60] + "\n")
    fh.write("//\n")


# ---------------------------------------------------------------------------
# the Y2H/MS parquet + the candidate-network mapping
# ---------------------------------------------------------------------------


def write_parquet_and_mapping(out, parquet_path, families, n_rows, min_n_tested, rng):
    """A real slice of the screen, made *dense*, plus a mapping that resolves its
    genes offline.

    BUILD_CANDIDATE_NETWORK's only network call is gene -> UniProt -> Pfam; given
    `--mapping-in` it makes none. Writing the mapping here from the genes the
    slice actually contains is what makes the fixture self-sufficient, and
    assigning our families to those genes is what makes the candidate network
    non-empty (both families of a pair must appear in a 3did DDI to survive).

    Non-empty is not enough, and this is the trap the density rows below exist
    for. `SAMPLE_NEGATIVES_ILP` draws a split's `ilp_candidates` negatives from
    the candidate pairs whose **both** endpoints landed in that split, and
    ppi-splitting's partitions are family-exclusive -- so survival is
    *quadratic* in the split's share, and the sampler hard-errors ("need N
    negatives but only M candidate pairs are available") rather than degrading.
    A real slice pairs whatever genes the screen happened to test, which after
    the round-robin family assignment collapsed to 34 pairs over 35 families:
    the 0.1 val split kept a median of **zero** usable pairs for ~4 positives.

    So the real rows are kept for provenance and every unordered family pair
    (self-pairs included -- they are legitimate DDIs) is appended on top, using
    one representative gene per family. BUILD_CANDIDATE_NETWORK subtracts the
    3did positives itself, so the generator does not need to know them: the
    emitted network is exactly the complete non-positive pair set over the
    fixture families, which is what the smallest split has to be sized against.
    """
    import pyarrow as pa
    import pyarrow.parquet as pq

    pf = pq.ParquetFile(parquet_path)
    columns = ["gene_name_bait", "gene_name_prey", "n_tested"]
    table = None
    for batch in pf.iter_batches(batch_size=65536, columns=columns):
        mask = [v is not None and v >= min_n_tested for v in batch.column("n_tested").to_pylist()]
        rows = [i for i, keep in enumerate(mask) if keep][:n_rows]
        if not rows:
            continue
        table = batch.take(rows)
        break
    if table is None:
        raise SystemExit(f"no parquet row has n_tested >= {min_n_tested}")

    real = pa.Table.from_batches([table])
    genes = sorted(
        {g for g in real.column("gene_name_bait").to_pylist() if g}
        | {g for g in real.column("gene_name_prey").to_pylist() if g}
    )
    if len(genes) < len(families):
        raise SystemExit(
            f"parquet slice has {len(genes)} genes but {len(families)} families to cover; "
            "raise --parquet-rows"
        )
    gene_to_uniprot, uniprot_to_pfams = {}, {}
    for i, gene in enumerate(genes):
        acc = f"C{i:05d}"
        gene_to_uniprot[gene] = acc
        uniprot_to_pfams[acc] = [families[i % len(families)]]

    # One representative gene per family: the first in `genes` order that the
    # round-robin above assigned to it, so the pair rows and the mapping cannot
    # disagree.
    rep = {}
    for i, gene in enumerate(genes):
        rep.setdefault(families[i % len(families)], gene)

    baits, preys = [], []
    for i, fam_a in enumerate(families):
        for fam_b in families[i:]:
            baits.append(rep[fam_a])
            preys.append(rep[fam_b])
    dense = pa.table(
        {
            "gene_name_bait": pa.array(baits, type=real.schema.field("gene_name_bait").type),
            "gene_name_prey": pa.array(preys, type=real.schema.field("gene_name_prey").type),
            # Every density row must clear params.negative_ppi_min_n_tested.
            "n_tested": pa.array([min_n_tested] * len(baits), type=real.schema.field("n_tested").type),
        }
    )

    path = os.path.join(out, "negative_ppi.parquet")
    pq.write_table(pa.concat_tables([real, dense]), path)
    log(f"wrote {path} ({real.num_rows} real rows + {dense.num_rows} density rows "
        f"covering all {len(families)} x {len(families)} family pairs)")

    mapping = os.path.join(out, "gene_pfam_mapping.json")
    with open(mapping, "w") as fh:
        json.dump({"gene_to_uniprot": gene_to_uniprot, "uniprot_to_pfams": uniprot_to_pfams}, fh, indent=1)
        fh.write("\n")  # keeps pre-commit's end-of-file-fixer from rewriting the fixture
    log(f"wrote {mapping} ({len(genes)} genes)")


# ---------------------------------------------------------------------------
# the enrichment inputs
# ---------------------------------------------------------------------------


def write_enrich_fixtures(out, records, families, swissprot_rows):
    """The UniProt / STRING / GO side, derived from the Pfam records.

    No local copy of these sources exists and the real ones are unusable as
    fixtures, so they are synthesised -- but from the accessions and coordinates
    the Pfam records actually carry, which is what keeps them consistent with what
    ingest will put in `domain_protein_map`.

    Embeddings are not fixtures at all any more: every model runs over the *cut
    domain sequence* on a GPU (ESM3/ESMC from gated weights, ProtT5 downloaded),
    and the results are published as HDF5 files rather than stored in the
    database. That is why `-profile test` is a cluster run. The synthetic
    `prott5.h5` this used to write went with the retired per-residue protein
    embeddings.
    """
    # accession -> (entry_name, [(start, end, domain_sequence)])
    proteins = {}
    for family in families:
        for header, seq in records[family]:
            m = PFAM_HEADER_RE.match(header[1:])
            entry, acc = m["entry"], m["acc"]
            start, end = int(m["start"]), int(m["end"])
            name, spans = proteins.setdefault(acc, (entry, []))
            spans.append((start, end, seq))

    # A parent sequence long enough to contain every one of its domains, with the
    # real domain sequences spliced in at their real coordinates -- so start/end
    # in domain_protein_map actually index the stored protein sequence.
    sequences = {}
    for acc, (_entry, spans) in proteins.items():
        length = max(end for _s, end, _q in spans)
        chars = ["A"] * length
        for start, end, domain_seq in spans:
            window = domain_seq[: end - start + 1].ljust(end - start + 1, "A")
            chars[start - 1 : end] = list(window)
        sequences[acc] = "".join(chars)

    # One flat file, three consumers: PARSE_SWISSPROT carves the protein FASTA,
    # the GO-term TSV and the accession->Pfam TSV back out of it. The two protein
    # sets stay as separate as they were when those were three fixture files --
    # Pfam-derived accessions carry GO and a sequence, HIPPIE-derived ones carry
    # the Pfam xref the single-domain step looks up -- so the sparser-than-real
    # fixture drives exactly the behaviour it drove before.
    dat = os.path.join(out, "uniprot_universe.dat.gz")
    hippie_side = {acc: (entry_name, gene, family)
                   for acc, entry_name, gene, family in swissprot_rows}
    with gzip.open(dat, "wt") as fh:
        for acc in sorted(set(proteins) | set(hippie_side)):
            entry_name, _spans = proteins.get(acc, (None, None))
            hippie = hippie_side.get(acc)
            if hippie:
                entry_name = entry_name or hippie[0]
            write_dat_entry(
                fh,
                acc=acc,
                entry_name=entry_name,
                gene=hippie[1] if hippie else entry_name.split("_")[0],
                seq=sequences.get(acc, "MSSRSVSRSR"),
                pfams=[hippie[2]] if hippie else [],
                go=["GO:0005515", "GO:0005829"] if acc in proteins else [],
            )
    log(f"wrote {dat} ({len(proteins)} Pfam-derived + {len(hippie_side)} HIPPIE-derived entries)")

    pfam2go = os.path.join(out, "pfam2go.txt")
    with open(pfam2go, "w") as fh:
        fh.write("!date: synthetic fixture\n")
        for family in families:
            fh.write(f"Pfam:{family} {family}_name > GO:protein binding ; GO:0005515\n")
    log(f"wrote {pfam2go} ({len(families)} families)")

    # STRING ids and the UniProt -> STRING mapping INSERT_PPI joins them through.
    ensp = {acc: f"9606.ENSP{i:011d}" for i, acc in enumerate(sorted(proteins))}
    mapping = os.path.join(out, "uniprot_id_mapping.dat.gz")
    with gzip.open(mapping, "wt") as fh:
        for acc, sid in sorted(ensp.items()):
            fh.write(f"{acc}\tSTRING\t{sid}\n")
            fh.write(f"{acc}\tGene_Name\t{proteins[acc][0].split('_')[0]}\n")
    log(f"wrote {mapping}")

    string = os.path.join(out, "string.txt.gz")
    ordered = [ensp[acc] for acc in sorted(proteins)]
    with gzip.open(string, "wt") as fh:
        fh.write("protein1 protein2 combined_score\n")
        for a, b in zip(ordered, ordered[1:]):
            fh.write(f"{a} {b} 900\n")
    log(f"wrote {string} ({max(len(ordered) - 1, 0)} edges)")


# ---------------------------------------------------------------------------


def choose_families(args, ddi_pfam_pairs, ppidm_pairs, negatome_pairs, fasta_path):
    degree = Counter()
    for a, b in ddi_pfam_pairs:
        degree[a] += 1
        degree[b] += 1
    external = Counter()
    for a, b in list(ppidm_pairs.keys()) + negatome_pairs:
        external[a] += 1
        external[b] += 1

    shortlist = sorted(
        degree,
        key=lambda f: (-min(degree[f], 20), -(1 if external[f] else 0), f),
    )[: args.shortlist]
    log(f"shortlisted {len(shortlist)} families for the Pfam-A.fasta pass")

    records, seen = scan_pfam_fasta(fasta_path, set(shortlist), args.max_records_per_family)
    eligible = [f for f in shortlist if seen[f] >= args.min_human_instances]
    eligible_set = set(eligible)
    log(f"{len(eligible)} of {len(shortlist)} shortlisted families have "
        f">= {args.min_human_instances} human instances")

    # Take the external pairs first. Both families of a pair have to be in the
    # fixture for the pair to survive PRUNE_UNREPRESENTED_DDIS, so picking
    # families one at a time by 3did degree tends to produce a fixture with no
    # external DDIs at all -- and then the external test set, the merge rules and
    # the conflict report are all untested. Negatome leads because it is the only
    # source of external *negatives*, and it is the smaller of the two.
    families = []
    for a, b in negatome_pairs + list(ppidm_pairs.keys()):
        if len(families) + 2 > args.families:
            break
        if a in eligible_set and b in eligible_set and a not in families and b not in families:
            families.extend([a, b] if a != b else [a])
    n_external = len(families)
    for family in eligible:
        if len(families) >= args.families:
            break
        if family not in families:
            families.append(family)
    if len(families) < 2:
        raise SystemExit("fewer than two families survived; lower --min-human-instances")
    log(f"families: {' '.join(families)} ({n_external} chosen for external-source coverage)")
    return families, records


def main():
    args = parse_args()
    rng = random.Random(args.seed)
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    out = args.out or os.path.join(repo, "tests", "data")
    os.makedirs(out, exist_ok=True)

    src = lambda name: os.path.join(args.sources, name)  # noqa: E731

    preamble, blocks, ddi_name_pairs, domain_rows = read_3did(src("3did.sql"))
    # Domain: (Name, oDB, ..., Pfam_id, ...) -- insert_3did joins DDI1.domain1/2
    # against Domain.Name and reads Pfam_id, so the fixture needs both columns.
    name_to_pfam = {}
    for row in domain_rows:
        acc = next((pfam_of(f) for f in row if pfam_of(f)), None)
        if acc:
            name_to_pfam[row[0]] = acc
    ddi_pfam_pairs = [
        (name_to_pfam[a], name_to_pfam[b])
        for a, b in ddi_name_pairs
        if a in name_to_pfam and b in name_to_pfam
    ]
    log(f"3did: {len(ddi_pfam_pairs)} DDIs resolve to Pfam accessions")

    ppidm_header, ppidm_pairs = read_ppidm(src("predicted_ddi_ppi_ppidm.tsv"))
    negatome_pairs = read_negatome(src("negatome_combined_pfam.txt"))

    families, records = choose_families(
        args, ddi_pfam_pairs, ppidm_pairs, negatome_pairs, src("Pfam-A.fasta.gz")
    )
    keep = set(families)

    # 3did rows over the kept families only.
    keep_domain = [row for row in domain_rows if name_to_pfam.get(row[0]) in keep]
    keep_names = {row[0] for row in keep_domain}
    keep_ddi_names = [(a, b) for a, b in ddi_name_pairs if a in keep_names and b in keep_names]
    keep_ddi_rows = [[a, "Pfam", b, "Pfam", "NULL", "a", "NULL", "NULL", "NULL"] for a, b in keep_ddi_names]
    log(f"3did fixture: {len(keep_ddi_rows)} DDIs over {len(keep_domain)} domain rows")
    write_3did_fixture(out, preamble, blocks, keep_ddi_rows, keep_domain)

    # External sources: pairs whose *both* families are kept, so they survive to
    # the external test set instead of being pruned for having no instances.
    ppidm_keep = [
        (a, b, cls)
        for (a, b), classes in sorted(ppidm_pairs.items())
        if a in keep and b in keep
        for cls in sorted(classes)
    ]
    with open(os.path.join(out, "ppidm.tsv"), "w") as fh:
        fh.write("domain_1\tdomain_2\tclass\n")
        for a, b, cls in ppidm_keep:
            fh.write(f"{a}\t{b}\t{cls}\n")
    log(f"wrote ppidm.tsv ({len(ppidm_keep)} rows)")

    negatome_keep = [(a, b) for a, b in negatome_pairs if a in keep and b in keep]
    with open(os.path.join(out, "negatome.txt"), "w") as fh:
        for a, b in negatome_keep:
            fh.write(f"{a}\t{b}\n")
    log(f"wrote negatome.txt ({len(negatome_keep)} pairs)")
    if not ppidm_keep and not negatome_keep:
        log("WARNING: no external pair survived, so the external test set will be empty")

    pfam_dir = write_pfam_fixtures(out, records, src("Pfam-A.clans.tsv"), families, args.pfam_release)
    mnemonics, accessions = set(), set()
    for family in families:
        for header, _seq in records[family]:
            m = PFAM_HEADER_RE.match(header[1:])
            mnemonics.add(m["entry"].rsplit("_", 1)[-1])
            accessions.add(m["acc"])
    write_interpro_cache(pfam_dir, args.pfam_release, mnemonics, accessions)

    hippie_rows, swissprot_rows = build_hippie_and_swissprot(
        src("HIPPIE-current.txt"), families, args.hippie_pairs, 0.63, rng
    )
    write_hippie(out, hippie_rows)

    write_parquet_and_mapping(
        out, src("negative_data_y2h_ms.pq"), families, args.parquet_rows, 5, rng
    )

    write_enrich_fixtures(out, records, families, swissprot_rows)

    # The three entries that make strata 1-3 of ppi-splitting's tier ladder
    # fillable: everything written above is human and Reviewed, so a test of
    # `--instance_tier all_species_reviewed` or `human_any_review_status` would
    # otherwise pass vacuously. Appended last, and to the Pfam-A.regions fixture
    # too -- which this script still does not generate, see the note in
    # tests/data/README.md. The default `human_reviewed` run ignores all three.
    for message in universe_tier_fixtures.apply(out):
        log(message)

    # Placeholders from the old `-stub` fixture set, superseded by the files
    # above, plus prott5.h5 from the retired per-residue protein embeddings.
    stale_prott5 = os.path.join(out, "prott5.h5")
    if os.path.exists(stale_prott5):
        os.remove(stale_prott5)
        log("removed prott5.h5 (per-residue protein embeddings are retired)")
    for stale in ("uniprot_go_terms.tsv", "negatome_combined_pfam.txt"):
        path = os.path.join(out, stale)
        if os.path.exists(path) and os.path.getsize(path) == 0:
            os.remove(path)
            log(f"removed empty placeholder {stale}")

    # Superseded by uniprot_universe.dat.gz, which PARSE_SWISSPROT splits into all three.
    for stale in ("swissprot_pfam.tsv", "uniprot_go_terms.tsv.gz",
                  "uniprot_sequences.fasta.gz"):
        path = os.path.join(out, stale)
        if os.path.exists(path):
            os.remove(path)
            log(f"removed {stale} (now carved out of uniprot_universe.dat.gz)")

    log("done")


if __name__ == "__main__":
    main()
