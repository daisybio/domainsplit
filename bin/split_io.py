#!/usr/bin/env python3
"""Readers for ppi-splitting's split CSVs, and the naming rules at that boundary.

ppi-splitting writes two levels of split CSV per (dataset row, negative set,
split):

* **family level** -- ``protein1,protein2,label`` plus whatever extra columns
  rode along from the positives file.  In DDI mode ``protein1``/``protein2`` are
  Pfam accessions; rows with ``label=0`` are the sampled negatives.
* **instance level** (``*_instances.csv``) -- ``protein1,protein2,label,
  family1,family2``, where ``protein1``/``protein2`` are *instance ids*
  (``{family}_{uniprot}_{start}_{end}``, the same key the ESM H5 files use) and
  ``family1``/``family2`` are the families they instantiate.

Neither file carries the split's identity: ppi-splitting puts it in the file
name, and the name only gains a negative-set suffix when a row asks for more
than one set.  Parsing names here would couple ingest to that suffixing rule, so
every ingest script instead takes explicit ``--split method:split:path`` triples
built by the calling subworkflow -- which is also where the negative-set ->
method-directory mapping (``ilp`` -> ``minimal_leakage``, ``ilp_candidates`` ->
``minimal_leakage_hcni``, ...) lives.

``val`` is renamed to ``validation`` here, at the ingest boundary; ppi-splitting
keeps ``val`` internally.
"""

import csv

#: ppi-splitting's split label -> the name domainsplit publishes it under.
SPLIT_RENAME = {"val": "validation"}

FAMILY_COLUMNS = ["protein1", "protein2", "label"]
INSTANCE_COLUMNS = ["protein1", "protein2", "label", "family1", "family2"]


def normalise_split(name):
    """``val`` -> ``validation``; every other label passes through unchanged."""
    return SPLIT_RENAME.get(name, name)


def parse_split_spec(spec):
    """``"minimal_leakage:val:/path/to/val.csv"`` -> ``("minimal_leakage", "validation", path)``.

    The path may itself contain colons, so only the first two are separators.
    """
    parts = spec.split(":", 2)
    if len(parts) != 3 or not all(parts):
        raise ValueError(f"--split expects method:split:path, got {spec!r}")
    method, split, path = parts
    return method, normalise_split(split), path


def _read_rows(path, required):
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = [c for c in required if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path} is missing column(s): {', '.join(missing)}")
        for row in reader:
            yield row


def read_family_csv(path):
    """Yield ``(pfam_a, pfam_b, label)`` from a family-level split CSV.

    An empty split still carries its header, so a zero-row file is valid input
    and yields nothing.
    """
    for row in _read_rows(path, FAMILY_COLUMNS):
        a, b = row["protein1"].strip(), row["protein2"].strip()
        if a and b:
            yield a, b, int(row["label"])


def read_instance_csv(path):
    """Yield ``(instance_a, instance_b, label, family_a, family_b)`` from a ``*_instances.csv``."""
    for row in _read_rows(path, INSTANCE_COLUMNS):
        ia, ib = row["protein1"].strip(), row["protein2"].strip()
        fa, fb = row["family1"].strip(), row["family2"].strip()
        if ia and ib and fa and fb:
            yield ia, ib, int(row["label"]), fa, fb


def read_instances_tsv(path):
    """Yield row dicts from ppi-splitting's ``instances.tsv``.

    Columns: ``instance_id, family, clan, protein_id, start, end, taxon_id,
    source_db`` -- the one table every DDI-mode consumer in ppi-splitting reads,
    kept whole here so the two pipelines cannot end up with two views of it.

    ``source_db`` is required, not optional: it carries the per-instance review
    status (``reviewed``/``unreviewed``) that ``ingest_instances.py`` stores as
    ``protein.reviewed``. Leaving it off the required list is how the flag went
    missing before -- the column existed upstream the whole time and this reader
    simply never asked for it.
    """
    columns = ["instance_id", "family", "clan", "protein_id", "start", "end", "taxon_id",
               "source_db"]
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        missing = [c for c in columns if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{path} is missing instances.tsv column(s): {', '.join(missing)}")
        for row in reader:
            yield row


def read_fasta(path):
    """``{header: sequence}`` for a plain FASTA.

    In DDI mode ppi-splitting keys ``sequences.fasta`` by *instance id* and the
    sequence is the domain's, not the parent protein's -- which is exactly what
    ``domain_protein_map.domain_sequence`` wants.
    """
    seqs, key, chunks = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if key is not None:
                    seqs[key] = "".join(chunks)
                key, chunks = line[1:].split()[0], []
            elif key is not None:
                chunks.append(line)
    if key is not None:
        seqs[key] = "".join(chunks)
    return seqs
