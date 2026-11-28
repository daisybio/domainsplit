#!/usr/bin/env python3
"""Write the computed-PPIDM fixtures under tests/data/ppidm/ (deterministic).

    tests/bin/make_ppidm_fixtures.py

Synthetic, not cut from the real sources: a few hundred PPIs over the 85
fixture families, shaped so that `-profile test` yields both label classes.

  * uniprot_ppidm.dat.gz   reviewed human entries X1P###, 1-2 `DR Pfam` each,
                           `DR STRING; 9606.FIXP#####`. Listed after
                           uniprot_universe.dat.gz in `ppidm_uniprot_dat`:
                           only 69 universe entries carry `DR Pfam`, one family
                           each, too few for a family to have the several
                           proteins a significant co-occurrence needs.
  * biogrid.mitab.txt.gz   BioGRID MITAB 2.5 (15 columns)
  * intact.mitab.txt.gz    IntAct MITAB 2.7 (42 columns; IntAct/MINT/DIP/UniProt)
  * hippie.txt.gz          every row of tests/data/hippie.tsv + planted rows
  * string/9606.protein.links.detailed.txt.gz
  * pdb_chain_uniprot.tsv.gz, 3did_flat.gz   chain pairs for the SIFTS source

Planted, among families with human instances and no 3did pair:
  * ENRICHED: two hetero pairs, every protein of one family bound to most
    proteins of the other in 7 of the 8 sources -> PPIDM Gold.
  * HUBS: four families with many high-degree proteins that never touch each
    other in any source -> every pair of them is powered and depleted in the
    six large sources -> INVERSE_PPIDM Gold. One protein carries two hubs, so
    a same-protein co-occurrence (which K must ignore) sits on a planted
    negative. Trap rows (IntAct Negative=true, BioGRID genetic, cross-taxon)
    also join hub proteins: an extractor that keeps them breaks the negatives.
  * gold-pair and random background PPIs, so 3did pairs get scored and the
    weight search / threshold scan have something to separate.
"""

import gzip
import random
import re
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "tests" / "data"
OUT = DATA / "ppidm"
SEED = 42
SOURCES = ["intact", "mint", "dip", "biogrid", "hippie", "string_exp", "string_rest", "sifts"]
ENRICHED_SOURCES = ["intact", "mint", "dip", "biogrid", "hippie", "string_exp", "sifts"]
HUB_SOURCES = ["intact", "mint", "biogrid", "hippie", "string_exp", "string_rest"]
BACKGROUND = {"intact": 70, "mint": 50, "dip": 15, "biogrid": 70, "hippie": 70,
              "string_exp": 60, "string_rest": 60, "sifts": 10}
HUB_PROTEINS, HUB_DEGREE = 6, 10
ENRICHED_PROTEINS = 4
GOLD_EDGES = 200
MULTI_DOMAIN = 10
METHODS = [("0018", "two hybrid"), ("0004", "affinity chromatography technology"),
           ("0096", "pull down"), ("0019", "coimmunoprecipitation")]
INTACT_DB = {"intact": ("0469", "IntAct"), "mint": ("0471", "MINT"), "dip": ("0465", "DIP")}


def write_gz(path, lines):
    """gzip with a fixed header (no name, mtime 0), so a rerun is byte-identical."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as gz:
        gz.write("".join(lines).encode())


def read_3did():
    """(pfam -> 3did name, set of canonical 3did Pfam pairs) from the fixture dump."""
    sql = gzip.open(DATA / "3did.sql.gz", "rt").read()
    dom = sql[sql.index("INSERT INTO `Domain`"):]
    dom = dom[:dom.index(";\n")]
    name_to_pf = {n: pf for n, pf in re.findall(r"\('([^']*)','Pfam',[^,]*,[^,]*,'(PF\d{5})", dom)}
    ddi = sql[sql.index("INSERT INTO `DDI1`"):]
    ddi = ddi[:ddi.index(";\n")]
    pairs = {tuple(sorted((name_to_pf[a], name_to_pf[b])))
             for a, b in re.findall(r"\('([^']*)','Pfam','([^']*)','Pfam'", ddi)
             if a in name_to_pf and b in name_to_pf}
    return {pf: n for n, pf in name_to_pf.items()}, pairs


def human_instance_families():
    """Families with a region on a reviewed human universe protein (instance_tier human_reviewed)."""
    human, cur = set(), {}
    with gzip.open(DATA / "uniprot_universe.dat.gz", "rt", encoding="latin-1") as fh:
        for line in fh:
            if line.startswith("ID"):
                cur = {"rev": "Reviewed;" in line.split()}
            elif line.startswith("AC") and "ac" not in cur:
                cur["ac"] = line[5:].split(";")[0].strip()
            elif line.startswith("OX"):
                cur["human"] = "NCBI_TaxID=9606;" in line
            elif line.startswith("//") and cur["rev"] and cur["human"]:
                human.add(cur["ac"])
    fams = set()
    with gzip.open(DATA / "pfam" / "Pfam-A.regions.tsv.gz", "rt") as fh:
        head = next(fh).rstrip("\n").split("\t")
        ia, ip = head.index("pfamseq_acc"), head.index("pfamA_acc")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[ia] in human:
                fams.add(f[ip])
    return fams


def universe_string_ids():
    """(STRING ids of universe entries with DR Pfam, STRING ids of those without)."""
    with_pf, without, sids, pf = [], [], [], False
    with gzip.open(DATA / "uniprot_universe.dat.gz", "rt", encoding="latin-1") as fh:
        for line in fh:
            if line.startswith("DR   STRING;"):
                sids.append(line.split(";")[1].strip())
            elif line.startswith("DR   Pfam;"):
                pf = True
            elif line.startswith("//"):
                (with_pf if pf else without).extend(sids)
                sids, pf = [], False
    return sorted(with_pf), sorted(without)


def pick_planted(rng, fams, gold):
    """Four mutually non-3did hub families, then two non-3did pairs disjoint from them."""
    order = sorted(fams)
    rng.shuffle(order)
    hubs = []
    for f in order:
        if len(hubs) < 4 and all(tuple(sorted((f, h))) not in gold for h in hubs):
            hubs.append(f)
    rest = [f for f in order if f not in hubs]
    enriched = []
    for a, b in combinations(rest, 2):
        used = {x for p in enriched for x in p}
        if len(enriched) < 2 and a not in used and b not in used and tuple(sorted((a, b))) not in gold:
            enriched.append(tuple(sorted((a, b))))
    assert len(hubs) == 4 and len(enriched) == 2
    return hubs, enriched


class Network:
    def __init__(self, rng, fams, hubs, enriched, gold):
        self.rng, self.hubs = rng, set(hubs)
        self.pfams = []                     # protein index -> sorted family list
        self.by_fam = {f: [] for f in fams}
        per_fam = {f: 2 for f in fams}
        per_fam.update({f: HUB_PROTEINS for f in hubs})
        per_fam.update({f: ENRICHED_PROTEINS for p in enriched for f in p})
        for f in sorted(fams):
            for _ in range(per_fam[f]):
                self.add([f])
        self.add(hubs[:2])                  # same-protein co-occurrence of a planted negative
        # A second domain only on plain proteins, so the planted families stay clean.
        plain = sorted(set(fams) - self.hubs - {f for p in enriched for f in p})
        candidates = [i for i, fs in enumerate(self.pfams) if set(fs) <= set(plain)]
        for i in rng.sample(candidates, MULTI_DOMAIN):
            second = rng.choice([f for f in plain if f not in self.pfams[i]])
            self.pfams[i] = sorted(self.pfams[i] + [second])
            self.by_fam[second].append(i)
        self.edges = {s: {} for s in SOURCES}   # source -> {(p, q): hi_conf}
        self.gold, self.enriched = gold, enriched

    def add(self, families):
        self.pfams.append(sorted(families))
        for f in families:
            self.by_fam[f].append(len(self.pfams) - 1)

    def hub_side(self, p):
        return bool(set(self.pfams[p]) & self.hubs)

    def link(self, source, p, q, hi):
        if p == q or (self.hub_side(p) and self.hub_side(q)):
            return False                    # hubs never touch each other: that is the negative
        key = (min(p, q), max(p, q))
        self.edges[source][key] = self.edges[source].get(key, False) or hi
        return True

    def build(self):
        rng, n = self.rng, len(self.pfams)
        for a, b in self.enriched:
            for s in ENRICHED_SOURCES:
                pairs = [(p, q) for p in self.by_fam[a] for q in self.by_fam[b]]
                for p, q in rng.sample(pairs, max(1, len(pairs) * 3 // 4)):
                    self.link(s, p, q, True)
        others = [p for p in range(n) if not self.hub_side(p)]
        for h in sorted(self.hubs):
            for p in self.by_fam[h]:
                for s in HUB_SOURCES:
                    for q in rng.sample(others, HUB_DEGREE):
                        self.link(s, p, q, rng.random() < 0.5)
        gold = sorted(g for g in self.gold if g[0] != g[1] and g[0] in self.by_fam and g[1] in self.by_fam)
        for a, b in rng.sample(gold, min(GOLD_EDGES, len(gold))):
            p, q = rng.choice(self.by_fam[a]), rng.choice(self.by_fam[b])
            for s in rng.sample(SOURCES, rng.randint(3, 6)):
                self.link(s, p, q, rng.random() < 0.8)
        for s, k in BACKGROUND.items():
            done = 0
            while done < k:
                done += self.link(s, rng.randrange(n), rng.randrange(n), rng.random() < 0.4)


def acc(p):
    return f"X1P{p:03d}"


def string_id(p):
    return f"9606.FIXP{p:05d}"


def gene(p):
    return f"FIX{p:03d}"


def write_dat(net):
    lines = []
    for p, fams in enumerate(net.pfams):
        seq = "M" + "A" * 59
        lines += [f"ID   {gene(p)}_HUMAN             Reviewed;          60 AA.\n",
                  f"AC   {acc(p)};" + (" X1S000;\n" if p == 0 else "\n"),
                  f"DE   RecName: Full=Synthetic PPIDM fixture protein {acc(p)};\n",
                  f"GN   Name={gene(p)};\n",
                  "OS   Homo sapiens (Human).\n",
                  "OX   NCBI_TaxID=9606;\n",
                  f"DR   STRING; {string_id(p)}; -.\n"]
        lines += [f"DR   Pfam; {f}; fixture; 1.\n" for f in fams]
        lines += ["SQ   SEQUENCE   60 AA;  0 MW;  0000000000000000 CRC64;\n", f"     {seq}\n", "//\n"]
    # An unreviewed entry with Pfam: not in the PPIDM universe, its PPIs are unmapped rows.
    lines += ["ID   FIXTR_HUMAN             Unreviewed;        60 AA.\n", "AC   X1T000;\n",
              "OX   NCBI_TaxID=9606;\n", f"DR   Pfam; {net.pfams[0][0]}; fixture; 1.\n",
              "SQ   SEQUENCE   60 AA;  0 MW;  0000000000000000 CRC64;\n", "     " + "A" * 60 + "\n", "//\n"]
    write_gz(OUT / "uniprot_ppidm.dat.gz", lines)


def mi(code, name):
    return f'psi-mi:"MI:{code}"({name})'


def write_biogrid(net, rng, traps):
    header = ["#ID Interactor A", "ID Interactor B", "Alt IDs Interactor A", "Alt IDs Interactor B",
              "Aliases Interactor A", "Aliases Interactor B", "Interaction Detection Method",
              "Publication 1st Author", "Publication Identifiers", "Taxid Interactor A",
              "Taxid Interactor B", "Interaction Types", "Source Database",
              "Interaction Identifiers", "Confidence Values"]
    rows, iid = ["\t".join(header) + "\n"], [100000]

    def row(p, q, pmid, itype=("0915", "physical association"), tax_b=9606, alt_b=None):
        iid[0] += 1
        # protein 0 is written by its secondary accession: it must resolve to X1P000
        side = lambda x: (f"biogrid:{200000 + x}|entrez gene/locuslink:{gene(x)}|"
                          f"uniprot/swiss-prot:{'X1S000' if x == 0 else acc(x)}|refseq:NP_{900000 + x}")
        code, name = rng.choice(METHODS)
        rows.append("\t".join([
            f"entrez gene/locuslink:{500000 + p}", f"entrez gene/locuslink:{500000 + q}",
            side(p), alt_b if alt_b is not None else side(q),
            f"entrez gene/locuslink:{gene(p)}(gene name)", f"entrez gene/locuslink:{gene(q)}(gene name)",
            mi(code, name), "Fixture A (2020)", f"pubmed:{pmid}", "taxid:9606", f"taxid:{tax_b}",
            mi(*itype), mi("0463", "biogrid"), f"biogrid:{iid[0]}", "-"]) + "\n")

    for (p, q), hi in sorted(net.edges["biogrid"].items()):
        physical = [("0915", "physical association"), ("0407", "direct interaction"),
                    ("0403", "colocalization"), ("0914", "association")]
        for k in range(2 if hi else 1):     # hi_conf = >= 2 distinct PubMed ids
            row(p, q, 30000000 + 1000 * p + q + k, rng.choice(physical))
        if hi and rng.random() < 0.3:
            row(p, q, 30000000 + 1000 * p + q)  # same PubMed id again: still two distinct
    p, q = traps
    row(p, q, 31000001, ("2373", "synthetic lethality"))   # genetic: dropped
    row(p, q, 31000002, ("0794", "synthetic genetic interaction defined by inequality"))
    row(p, q, 31000003, tax_b=10090)                       # cross-taxon: dropped
    row(p, q, 31000004, alt_b=f"biogrid:1|entrez gene/locuslink:{gene(q)}")  # no swiss-prot: dropped
    row(p, q, 31000005, alt_b=f"uniprot/swiss-prot:{acc(q)}|uniprot/swiss-prot:{acc(p)}")  # ambiguous
    write_gz(OUT / "biogrid.mitab.txt.gz", rows)


def write_intact(net, rng, traps):
    header = ["#ID(s) interactor A", "ID(s) interactor B", "Alt. ID(s) interactor A",
              "Alt. ID(s) interactor B", "Alias(es) interactor A", "Alias(es) interactor B",
              "Interaction detection method(s)", "Publication 1st author(s)",
              "Publication Identifier(s)", "Taxid interactor A", "Taxid interactor B",
              "Interaction type(s)", "Source database(s)", "Interaction identifier(s)",
              "Confidence value(s)", "Expansion method(s)", "Biological role(s) interactor A",
              "Biological role(s) interactor B", "Experimental role(s) interactor A",
              "Experimental role(s) interactor B", "Type(s) interactor A", "Type(s) interactor B",
              "Xref(s) interactor A", "Xref(s) interactor B", "Interaction Xref(s)",
              "Annotation(s) interactor A", "Annotation(s) interactor B",
              "Interaction annotation(s)", "Host organism(s)", "Interaction parameter(s)",
              "Creation date", "Update date", "Checksum(s) interactor A", "Checksum(s) interactor B",
              "Interaction Checksum(s)", "Negative", "Feature(s) interactor A",
              "Feature(s) interactor B", "Stoichiometry(s) interactor A",
              "Stoichiometry(s) interactor B", "Identification method participant A",
              "Identification method participant B"]
    assert len(header) == 42
    rows, iid = ["\t".join(header) + "\n"], [0]
    human = "taxid:9606(human)|taxid:9606(Homo sapiens)"

    def row(id_a, id_b, db, score, negative="false", tax_b=human, type_b=mi("0326", "protein")):
        iid[0] += 1
        code, name = rng.choice(METHODS)
        f = [id_a, id_b, f"intact:EBI-{iid[0]}1", f"intact:EBI-{iid[0]}2", "-", "-", mi(code, name),
             "Fixture et al. (2021)", f"pubmed:{32000000 + iid[0]}|imex:IM-{iid[0]}-1", human, tax_b,
             mi("0915", "physical association"), mi(*db), f"intact:EBI-{iid[0]}0",
             f"author score:high|intact-miscore:{score:.2f}",
             mi("1060", "spoke expansion") if iid[0] % 2 else "-",
             mi("0499", "unspecified role"), mi("0499", "unspecified role"),
             mi("0496", "bait"), mi("0498", "prey"), mi("0326", "protein"), type_b,
             "-", "-", "-", "-", "-", "-", "taxid:-1(in vitro)", "-", "2021/01/01", "2026/01/01",
             "-", "-", "-", negative, "-", "-", "-", "-",
             mi("0396", "predetermined participant"), mi("0396", "predetermined participant")]
        assert len(f) == 42
        rows.append("\t".join(f) + "\n")

    for source in ("intact", "mint", "dip"):
        for (p, q), hi in sorted(net.edges[source].items()):
            # any database other than MINT/DIP counts as intact
            db = ("0486", "UniProt") if source == "intact" and rng.random() < 0.1 else INTACT_DB[source]
            # isoform / chain suffixes are stripped by the extractor
            suffix = rng.choice(["", "", "", "-2", "-PRO_0000012345"])
            score = rng.uniform(0.45, 0.9) if hi else rng.uniform(0.2, 0.4)  # hi: > 0.4
            row(f"uniprotkb:{acc(p)}{suffix}", f"uniprotkb:{acc(q)}", db, score)
    p, q = traps
    row(f"uniprotkb:{acc(p)}", f"uniprotkb:{acc(q)}", INTACT_DB["intact"], 0.9, negative="true")
    row(f"uniprotkb:{acc(p)}", f"uniprotkb:{acc(q)}", INTACT_DB["mint"], 0.9,
        tax_b="taxid:10090(mouse)|taxid:10090(Mus musculus)")
    row(f"uniprotkb:{acc(p)}", 'chebi:"CHEBI:15422"', INTACT_DB["intact"], 0.5,
        type_b=mi("0328", "small molecule"))
    row(f"uniprotkb:{acc(p)}", "uniprotkb:X1T000", INTACT_DB["intact"], 0.5)  # unreviewed: unmapped
    write_gz(OUT / "intact.mitab.txt.gz", rows)


def write_hippie(net, rng):
    lines = (DATA / "hippie.tsv").read_text().splitlines(keepends=True)
    assert lines[0].startswith("uniprot_accession_A\t")
    for (p, q), hi in sorted(net.edges["hippie"].items()):
        # 0.72 itself is not high-confidence (strict >)
        score = rng.choice([0.73, 0.7451, 0.82, 0.9]) if hi else rng.choice([0.3, 0.5, 0.63, 0.72])
        lines.append(f"{acc(p)}\t{gene(p)}_HUMAN\t{500000 + p}\t{acc(q)}\t{gene(q)}_HUMAN\t"
                     f"{500000 + q}\t{score:.4f}\texperiments:Two-hybrid;pmids:{33000000 + p};"
                     f"sources:biogrid\n")
    write_gz(OUT / "hippie.txt.gz", lines)


def write_string(net, rng):
    pairs = {}                              # (a, b) -> [exp, rest, hi]
    for name, i in (("string_exp", 0), ("string_rest", 1)):
        for k, hi in net.edges[name].items():
            e = pairs.setdefault((string_id(k[0]), string_id(k[1])), [False, False, False])
            e[i], e[2] = True, e[2] or hi
    with_pf, without = universe_string_ids()
    for i, sid in enumerate(with_pf + without[:3]):  # universe ids; Pfam-less ones are unmapped
        pairs[(sid, string_id(i))] = [True, True, False]
    rows = ["protein1 protein2 neighborhood fusion cooccurence coexpression experimental "
            "database textmining combined_score\n"]
    for (a, b), (exp, rest, hi) in sorted(pairs.items()):
        comb = rng.randint(901, 999) if hi else rng.randint(150, 900)  # 900 itself is low
        e = (comb if not rest else rng.randint(100, comb - 50)) if exp else 0
        t = comb - e if rest else 0
        for x, y in ((a, b), (b, a)):       # STRING lists both directions
            rows.append(f"{x} {y} 0 0 0 0 {e} 0 {t} {comb}\n")
    write_gz(OUT / "string" / "9606.protein.links.detailed.txt.gz", rows)


def write_sifts(net, names, traps):
    sifts = ["# 2026/11/01 - 12:00 | PDB: 45.26 | UniProt: 2026.04\n",
             "PDB\tCHAIN\tSP_PRIMARY\tRES_BEG\tRES_END\tPDB_BEG\tPDB_END\tSP_BEG\tSP_END\n"]
    flat = []
    for i, (p, q) in enumerate(sorted(net.edges["sifts"])):
        pdb = f"9x{i:02d}"
        for chain, x in (("A", p), ("B", q)):
            sifts.append(f"{pdb}\t{chain}\t{acc(x)}\t1\t60\t1\t60\t1\t60\n")
        fa, fb = net.pfams[p][0], net.pfams[q][0]
        flat += [f"#=ID\t{names.get(fa, fa)}\t{names.get(fb, fb)}\t ({fa}.1@Pfam {fb}.1@Pfam)\n",
                 # upper-case PDB id here, lower-case in SIFTS: the join normalises case
                 f"#=3D\t{pdb.upper() if i % 3 == 0 else pdb}\tA:1-60\tB:1-60\t12.0\t3.0\t1:1:1\n",
                 "//\n"]
    p, q = traps
    flat += ["#=ID\thubA\thubB\t (PF00000.1@Pfam PF00001.1@Pfam)\n",
             "#=3D\t9zz0\tA:1-60\tA:1-60\t12.0\t3.0\t1:1:1\n",      # same chain: dropped
             "#=3D\t9zz1\tA:1-60\tB:1-60\t12.0\t3.0\t1:1:1\n", "//\n"]  # chain B unmapped
    sifts.append(f"9zz0\tA\t{acc(p)}\t1\t60\t1\t60\t1\t60\n")
    sifts.append(f"9zz1\tA\t{acc(q)}\t1\t60\t1\t60\t1\t60\n")
    write_gz(OUT / "pdb_chain_uniprot.tsv.gz", sifts)
    write_gz(OUT / "3did_flat.gz", flat)


def main():
    rng = random.Random(SEED)
    names, gold = read_3did()
    fams = human_instance_families()
    assert len(fams) == 85, len(fams)
    hubs, enriched = pick_planted(rng, fams, gold)
    net = Network(rng, fams, hubs, enriched, gold)
    net.build()
    # traps join the two planted-negative hubs: kept by mistake, they would give K > 0
    traps = (net.by_fam[hubs[2]][0], net.by_fam[hubs[3]][0])
    write_dat(net)
    write_biogrid(net, rng, traps)
    write_intact(net, rng, traps)
    write_hippie(net, rng)
    write_string(net, rng)
    write_sifts(net, names, traps)
    print("hub families (every pair should be INVERSE_PPIDM Gold):", " ".join(hubs))
    print("enriched pairs (PPIDM Gold):", " ".join(f"{a}-{b}" for a, b in enriched))
    print("PPIs per source:", {s: len(e) for s, e in net.edges.items()})
    print("proteins:", len(net.pfams))


if __name__ == "__main__":
    main()
