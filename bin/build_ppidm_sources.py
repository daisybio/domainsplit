#!/usr/bin/env python3
"""Build the eight per-source PPI networks PPIDM's statistics run on.

Universe: every reviewed UniProt entry (all species) with >= 1 Pfam domain.
Writes into --outdir:

  protein_pfam.tsv.gz        acc pfams            (accession -> Pfam map)
  source_<name>.tsv.gz       acc_a acc_b hi_conf  (one per source, unique canonical pairs)
  ppidm_sources.tsv          source raw_rows mapped_pairs hi_conf_pairs unmapped_rows

Pairs are int64-coded (index of acc_a * M + index of acc_b over the sorted
universe), so canonical string order is canonical code order and the whole run
fits in a few GB despite tens of millions of STRING rows.
"""

import argparse
import gzip
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd

from parse_single_domain_ppi import detect_layout
from parse_swissprot_dat import TAXID_RE, review_flag

SOURCES = ["intact", "mint", "dip", "biogrid", "hippie", "string_exp", "string_rest", "sifts"]
INTACT_SPLIT = {"0471": "mint", "0465": "dip"}  # every other source_db -> intact


def log(msg):
    print(f"[build_ppidm_sources] {msg}", flush=True)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--sprot-dat", action="append", required=True)
    p.add_argument("--biogrid", required=True)
    p.add_argument("--intact", required=True)
    p.add_argument("--hippie", required=True)
    p.add_argument("--string-detailed", action="append", required=True)
    p.add_argument("--sifts-chain-uniprot", required=True)
    p.add_argument("--did-flat", required=True)
    p.add_argument("--hippie-min-score", type=float, required=True)
    p.add_argument("--string-min-score", type=float, required=True)
    p.add_argument("--intact-min-miscore", type=float, required=True)
    p.add_argument("--biogrid-min-pubs", type=int, required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def is_gzip(path):
    with open(path, "rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def open_text(path):
    return (gzip.open if is_gzip(path) else open)(path, "rt", encoding="latin-1")


def read_universe(paths):
    """Reviewed entries with Pfam, the secondary->primary map and the STRING map."""
    proteins = {}        # primary -> (taxid, sorted pfams, string ids)
    reviewed = set()     # every reviewed primary, Pfam or not
    sec_owner = {}       # secondary -> primary, or None when several entries list it
    for path in paths:
        log(f"reading {path}")
        accs, ox, pfams, strings, skip, keep = [], "", set(), [], False, False
        with open_text(path) as fh:
            for line in fh:
                tag = line[:2]
                if tag == "//":
                    if keep and accs and accs[0] not in reviewed:  # first file wins
                        primary = accs[0]
                        reviewed.add(primary)
                        for s in accs[1:]:
                            sec_owner[s] = primary if sec_owner.get(s, primary) == primary else None
                        if pfams:
                            proteins[primary] = (ox, sorted(pfams), strings)
                    accs, ox, pfams, strings, skip, keep = [], "", set(), [], False, False
                elif skip:
                    continue
                elif tag == "ID":
                    keep = review_flag(line, path)
                    skip = not keep
                elif tag == "AC":
                    accs.extend(t.strip() for t in line[5:].split(";") if t.strip())
                elif tag == "OX":
                    m = TAXID_RE.search(line)
                    ox = m.group(1) if m else ""
                elif tag == "DR":
                    body = line[5:]
                    if body.startswith("Pfam; "):
                        pfams.add(body[6:].split(";")[0].strip())
                    elif body.startswith("STRING; "):
                        strings.append(body[8:].split(";")[0].strip())
                elif tag == "SQ":
                    skip = True  # the sequence block runs to "//"
    return proteins, reviewed, sec_owner


class Universe:
    def __init__(self, proteins, reviewed, sec_owner):
        self.accs = np.array(sorted(proteins), dtype=object)
        self.M = len(self.accs)
        self.idx = {a: i for i, a in enumerate(self.accs)}
        self.taxid = np.array([proteins[a][0] for a in self.accs], dtype=object)
        # A secondary resolves only when one entry lists it and it is nobody's primary.
        self.lookup = dict(self.idx)
        n_sec = n_shared = 0
        for s, p in sec_owner.items():
            if p is None:
                n_shared += 1
            elif s not in reviewed and p in self.idx:
                self.lookup[s] = self.idx[p]
                n_sec += 1
        string_owner = {}
        for a in self.accs:
            for sid in proteins[a][2]:
                string_owner.setdefault(sid, set()).add(self.idx[a])
        self.string = {sid: next(iter(o)) for sid, o in string_owner.items() if len(o) == 1}
        self.n_string_ambiguous = len(string_owner) - len(self.string)
        log(f"universe: {self.M} reviewed proteins with Pfam; {n_sec} secondary accessions "
            f"resolve, {n_shared} shared by several entries do not; {len(self.string)} STRING ids, "
            f"{self.n_string_ambiguous} mapping to several accessions skipped")

    def codes(self, accs):
        """Accession strings (isoform/chain suffixes allowed) -> universe index, -1 if absent."""
        s = pd.Series(accs, dtype=object).str.split("-", n=1).str[0]
        return s.map(self.lookup).fillna(-1).to_numpy(dtype=np.int64)


def pair_keys(a, b, M):
    """Mask of rows with both sides mapped and distinct, and their canonical keys."""
    ok = (a >= 0) & (b >= 0) & (a != b)
    lo, hi = np.minimum(a[ok], b[ok]), np.maximum(a[ok], b[ok])
    return ok, lo * M + hi


def max_by_key(keys, values):
    s = pd.Series(values).groupby(keys).max()
    return s.index.to_numpy(dtype=np.int64), s.to_numpy()


class Writer:
    def __init__(self, uni, outdir):
        self.uni, self.outdir, self.rows = uni, outdir, []

    def write(self, name, keys, hi, raw, unmapped, note=""):
        keys = np.asarray(keys, dtype=np.int64)
        hi = np.asarray(hi, dtype=bool)
        df = pd.DataFrame({"acc_a": self.uni.accs[keys // self.uni.M],
                           "acc_b": self.uni.accs[keys % self.uni.M],
                           "hi_conf": hi.astype(int)})
        df.to_csv(os.path.join(self.outdir, f"source_{name}.tsv.gz"), sep="\t", index=False,
                  compression="gzip")
        self.rows.append((name, int(raw), len(keys), int(hi.sum()), int(unmapped)))
        log(f"{name}: raw_rows={raw} mapped_pairs={len(keys)} hi_conf_pairs={int(hi.sum())} "
            f"unmapped_rows={unmapped}{note}")


def mitab_sources(path, uni, w, kind, args):
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    a, b = uni.codes(df["acc_a"]), uni.codes(df["acc_b"])
    unmapped = (a < 0) | (b < 0)
    selfp = ~unmapped & (a == b)
    if kind == "biogrid":
        ok, keys = pair_keys(a, b, uni.M)
        pm = pd.DataFrame({"key": keys, "pmid": df["pmids"].to_numpy()[ok]})
        pm = pm.assign(pmid=pm["pmid"].str.split(";")).explode("pmid")
        pm = pm[pm["pmid"].fillna("") != ""].drop_duplicates()
        uniq = np.unique(keys)
        npubs = pm.groupby("key").size().reindex(uniq, fill_value=0).to_numpy()
        w.write("biogrid", uniq, npubs >= args.biogrid_min_pubs, len(df), unmapped.sum(),
                f" self_after_mapping={selfp.sum()}")
        return
    label = df["source_db"].map(INTACT_SPLIT).fillna("intact").to_numpy()
    score = pd.to_numeric(df["score"], errors="coerce").to_numpy()
    for name in ("intact", "mint", "dip"):
        sel = label == name
        ok, keys = pair_keys(a[sel], b[sel], uni.M)
        keys, best = max_by_key(keys, score[sel][ok])
        # NaN (no miscore) compares False: never high-confidence
        w.write(name, keys, best > args.intact_min_miscore, sel.sum(), unmapped[sel].sum(),
                f" self_after_mapping={selfp[sel].sum()}")


def hippie_source(path, uni, w, args):
    with open_text(path) as fh:
        ca, cb, cs = detect_layout(fh.readline())
        df = pd.read_csv(fh, sep="\t", header=None, usecols=[ca, cb, cs], dtype=str,
                         keep_default_na=False, quoting=3)
    a, b = uni.codes(df[ca]), uni.codes(df[cb])
    score = pd.to_numeric(df[cs], errors="coerce").to_numpy()
    ok, keys = pair_keys(a, b, uni.M)
    keys, best = max_by_key(keys, score[ok])
    w.write("hippie", keys, best > args.hippie_min_score, len(df), ((a < 0) | (b < 0)).sum())


def string_sources(paths, uni, w, args):
    parts = {"string_exp": [], "string_rest": []}
    raw, unmapped = Counter(), Counter()
    for path in paths:
        comp = "gzip" if is_gzip(path) else None
        cols = list(pd.read_csv(path, sep=" ", nrows=0, compression=comp).columns)
        exp = next((c for c in ("experimental", "experiments") if c in cols), None)
        if exp is None or not {"protein1", "protein2", "combined_score"} <= set(cols):
            raise SystemExit(f"ERROR: {path} is not a STRING protein.links.detailed file "
                             f"(columns: {cols})")
        per_file = {k: ([], []) for k in parts}
        for chunk in pd.read_csv(path, sep=" ", compression=comp, chunksize=2_000_000,
                                 usecols=["protein1", "protein2", exp, "combined_score"],
                                 dtype={"protein1": str, "protein2": str}):
            a = chunk["protein1"].map(uni.string).fillna(-1).to_numpy(dtype=np.int64)
            b = chunk["protein2"].map(uni.string).fillna(-1).to_numpy(dtype=np.int64)
            comb = chunk["combined_score"].to_numpy()
            e = chunk[exp].to_numpy()
            mapped = (a >= 0) & (b >= 0)
            for name, sel in (("string_exp", e > 0), ("string_rest", comb - e > 0)):
                raw[name] += int(sel.sum())
                unmapped[name] += int((sel & ~mapped).sum())
                ok, keys = pair_keys(np.where(sel, a, -1), b, uni.M)
                per_file[name][0].append(keys)
                per_file[name][1].append(comb[ok])
        for name, (k, v) in per_file.items():  # both directions sit in one file: dedupe early
            parts[name].append(max_by_key(np.concatenate(k), np.concatenate(v)))
        log(f"read {path}")
    for name, kv in parts.items():
        keys, best = max_by_key(np.concatenate([k for k, _ in kv]), np.concatenate([v for _, v in kv]))
        w.write(name, keys, best > args.string_min_score, raw[name], unmapped[name])


def sifts_source(sifts_path, did_path, uni, w):
    with open_text(sifts_path) as fh:
        first = fh.readline()
        header = fh.readline() if first.startswith("#") else first
        df = pd.read_csv(fh, sep="\t", header=None, names=header.rstrip("\n").split("\t"),
                         usecols=["PDB", "CHAIN", "SP_PRIMARY"], dtype=str, keep_default_na=False)
    codes = uni.codes(df["SP_PRIMARY"])
    chain_accs = {}
    for pdb, chain, c in zip(df["PDB"].str.lower(), df["CHAIN"], codes):
        if c >= 0:
            chain_accs.setdefault((pdb, chain), set()).add(int(c))

    chain_pairs = set()
    with open_text(did_path) as fh:
        for line in fh:
            if line.startswith("#=3D"):
                f = line.split("\t")
                ca, cb = f[2].split(":", 1)[0], f[3].split(":", 1)[0]
                if ca != cb:
                    chain_pairs.add((f[1].lower(), *sorted((ca, cb))))
    keys, unmapped, cross = set(), 0, 0
    for pdb, ca, cb in chain_pairs:
        sa, sb = chain_accs.get((pdb, ca)), chain_accs.get((pdb, cb))
        if not sa or not sb:
            unmapped += 1
            continue
        for x in sa:
            for y in sb:
                if x == y:
                    continue
                if uni.taxid[x] != uni.taxid[y]:
                    cross += 1
                    continue
                keys.add(min(x, y) * uni.M + max(x, y))
    keys = np.array(sorted(keys), dtype=np.int64)
    w.write("sifts", keys, np.ones(len(keys), dtype=bool), len(chain_pairs), unmapped,
            f" cross_taxon_accession_pairs={cross}")


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    proteins, reviewed, sec_owner = read_universe(args.sprot_dat)
    if not proteins:
        raise SystemExit("ERROR: no reviewed entry with a Pfam cross-reference in --sprot-dat")
    with gzip.open(os.path.join(args.outdir, "protein_pfam.tsv.gz"), "wt", compresslevel=6) as out:
        out.write("acc\tpfams\n")
        for a in sorted(proteins):
            out.write(f"{a}\t{';'.join(proteins[a][1])}\n")
    uni = Universe(proteins, reviewed, sec_owner)
    del proteins, reviewed, sec_owner

    w = Writer(uni, args.outdir)
    mitab_sources(args.intact, uni, w, "intact", args)
    mitab_sources(args.biogrid, uni, w, "biogrid", args)
    hippie_source(args.hippie, uni, w, args)
    string_sources(args.string_detailed, uni, w, args)
    sifts_source(args.sifts_chain_uniprot, args.did_flat, uni, w)

    order = {n: i for i, n in enumerate(SOURCES)}
    with open(os.path.join(args.outdir, "ppidm_sources.tsv"), "w") as f:
        f.write("source\traw_rows\tmapped_pairs\thi_conf_pairs\tunmapped_rows\n")
        for row in sorted(w.rows, key=lambda r: order[r[0]]):
            f.write("\t".join(map(str, row)) + "\n")

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    numpy: {np.__version__}\n")
        f.write(f"    pandas: {pd.__version__}\n")


if __name__ == "__main__":
    main()
