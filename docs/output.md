# daisybio/domainsplit: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and produces:

- [Master database](#master-database) - the single enriched SQLite database
- [Split databases](#split-databases) - one database per (splitting method, split)
- [Domain embeddings](#domain-embeddings) - one HDF5 per embedding model
- [Reports](#reports) - what left the database, or failed to enter it
- [Splitting inputs and diagnostics](#splitting-inputs-and-diagnostics)
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Master database

<details markdown="1">
<summary>Output files</summary>

- `domainsplit.sqlite3`
  - Every DDI that survived pruning, its domains and domain instances, protein sequences, GO
    annotations, the STRING PPI network, and `ddi_split_membership` — which records for every
    (method, split) exactly which DDIs and which concrete domain-instance pairs belong to it.

</details>

The master carries **no BLOBs**: embeddings are published separately (see below). `start_pos` and
`end_pos` on `domain_protein_map` are `INTEGER`, and `instance_id` is ppi-splitting's
`{pfam_id}_{uniprot_id}_{start}_{end}` — the same key the embedding files use.

### Split databases

<details markdown="1">
<summary>Output files</summary>

- `databases/<method>/<split>.sqlite3`
  - One database per (method, split), each a self-contained subset of the master holding only its
    own split's DDIs and everything they reach.
  - Only `ddi_split_membership` is the split assignment. `domain_protein_map`, `protein` and PPI
    rows of a kept family can include proteins ppi-splitting placed in another split, so do not
    form new evaluation pairs from them.

</details>

With the default `--ppi_splitting_multi_negset true` there are **7 method directories and 26
databases**:

| method                 | splits                                                         |
| ---------------------- | -------------------------------------------------------------- |
| `random`               | `train`, `validation`, `test_uniform_1_1`, `test_uniform_1_10` |
| `random_ilp`           | `train`, `validation`, `test_matched_1_1`, `test_uniform_1_10` |
| `random_hrni`          | `train`, `validation`, `test_matched_1_1`, `test_uniform_1_10` |
| `minimal_leakage`      | `train`, `validation`, `test_matched_1_1`, `test_uniform_1_10` |
| `minimal_leakage_hrni` | `train`, `validation`, `test_matched_1_1`, `test_uniform_1_10` |
| `external_test`        | `train`, `validation`, `test`                                  |
| `external_test_hrni`   | `train`, `validation`, `test`                                  |

`random` is the naive baseline (random split, uniform negatives). `minimal_leakage` is the main
product (ILP split, ILP negatives). `random_ilp` / `random_hrni` put the ILP negative sets on the
random split, so HRNI can be compared without the family-disjoint split as a confounder. The `_hrni` directories are the same positive split with
**h**igh-**r**eplication **n**on-**i**nteractions: negatives drawn from the candidate network of
Pfam pairs a Y2H/MS screen tested without finding an interaction. `external_test`'s `test` split is
not an internal held-out slice at all — it is the external DDIs (PPIDM / single-domain-PPI /
Negatome), unseen by construction.

The two internal test sets are named for the sampler that drew their negatives, because they
are not one test set at two class ratios. `test_uniform_1_10` is always drawn uniformly at 1:10.
The 1:1 set is drawn by the negative set's own sampler: uniform for `random`
(`test_uniform_1_1`), the bias-matched ILP for the `ilp`/`ilp_candidates` sets (`test_matched_1_1`). Comparing a
predictor's ROC-AUC between `test_matched_1_1` and `test_uniform_1_10` therefore compares two
negative distributions, not two ratios. A sampled negative that 3did already holds as a
positive is dropped (`reports/unresolved_split_rows.tsv`, `reason=label_mismatch`), so a 1:1 set
can be a few negatives short of exactly 1:1.

Set `--ppi_splitting_multi_negset false` and this drops to 3 directories / 11 databases.

### Domain embeddings

<details markdown="1">
<summary>Output files</summary>

- `embeddings/esm3_domain_embeddings.h5`
- `embeddings/esmc_domain_embeddings.h5`
- `embeddings/prott5_domain_embeddings.h5`

</details>

One mean-pooled `float16` vector per **domain instance** per model, keyed

```
h5[pfam_id][instance_id]
```

where `pfam_id` is the Pfam accession (`PF00069`) and `instance_id` is
`domain_protein_map.instance_id` — precisely `h5[pfam_id][COALESCE(instance_id, 'r' || rowid)]` over
the split database's `domain_protein_map`, joined to `domain` for the accession.

Only the **cut domain sequence** is embedded, never the parent protein. A per-residue protein
embedding sliced to a domain's coordinates still carries protein context, and protein context
correlates with the interaction partner — which would smuggle protein identity into a split that
was partitioned on families.

Root attributes on each file: `model`, `pooling`, `dim`, `dtype`, `key_layout`, `n_domains`,
`n_instances`.

> [!NOTE]
> The outer key used to be `domain.id`, a **surrogate integer**, which made a file valid across every
> split database of the same run and silently wrong across runs — the failure mode being a successful
> lookup of the wrong domain's vectors. `domain` is `UNIQUE(pfam_id)`, so the accession is the same
> information with run-independent meaning. There is deliberately no run identifier on the file: with
> accession keys a lookup either finds the family or does not, so pairing an embedding file with the
> databases it belongs to is the caller's business.

Which models run is `--embedding_models_domainbench` (default `esm3,esmc,prott5`); set it empty
(`--embedding_models_domainbench=`, or a params file) to skip embedding entirely. ESM3 and ESMC need `HF_TOKEN`; ProtT5 does not.

### Structures

<details markdown="1">
<summary>Output files</summary>

- `structures/structures.h5`
- `structures/structure_coverage.tsv`
- `structures/af3_to_model/<method>.tsv`, `structures/af3_to_model/all_methods.tsv`

</details>

Written only for the methods in `--structure_data_methods`.

`structures.h5` holds each instance pair's sliced two-chain PDB as gzip bytes at
`{pfam_id_a}_{pfam_id_b}/{instance_id_a}__{instance_id_b}`. Each slice dataset carries these attrs:

- `iptm`, `ptm`, `ranking_score`: float, NaN when absent.
- `af_source`: `inferred` for own AF3 runs, `downloaded` for AlphaFold DB complexes.
- `chain_pair_iptm`: 2×2, AF3 only, in the row's orientation (index 0 is the row's protein A).

AlphaFold DB models have no confidence summary. Their `ptm` and `ranking_score` are always NaN, and
their `iptm` is NaN unless the PDB header carries one. pLDDT is in the PDB B-factor column.

The PAE submatrix over the two sliced domain windows is at
`pae/{pfam_id_a}_{pfam_id_b}/{instance_id_a}__{instance_id_b}`:

- float32 and square, with rows and columns ordered [domain A residues, domain B residues].
- Rows are aligned residues and columns are scored residues, as in AlphaFold.
- Attrs `residues_a` and `residues_b` give the UniProt residue numbers.
- The dataset is absent when the PAE file is missing.

A confidence file that is missing is not fatal: it leaves NaN attrs and no PAE. One that is present
but cannot be parsed, or does not fit its model, stops the run.

`af_metadata/` is the AF model index this run used, rebuilt by `BUILD_AF_METADATA` from
`--af_downloaded_manifest`, `--af_inferred_dirs` and/or `--af_metadata` on every run:

- `af_metadata.csv`: one model per unordered protein pair, `uniprot_id_a` being chain A of the
  file, plus a stat fingerprint (`*_size`, `*_mtime_ns`; -1 = file absent) of the model, its PAE
  and its summary file. A changed fingerprint re-runs the slicing.
- `af_inventory.tsv`: every candidate with `status` `selected`, `shadowed` (an inferred model
  of a pair that has a downloaded one), `model_missing` or `exact_duplicate`, and its `origin`.
- `af_metadata_report.tsv`: those counts per source.

`structure_coverage.tsv` has one row per (`method`, `split`, `label`), with columns
`instance_pairs`, `with_structure`, `fraction`, `with_pae`, `with_iptm`, `same_protein_pairs` and
`same_protein_with_structure`. A same-protein pair has both instances on one protein; its structure
comes from a homodimer model, so the slice is an inter-chain contact even where the real one may be
intra-chain. Positives and negatives
are covered at very different rates, which makes structure features a label leak, so check this
file before using them.

`af3_to_model/` lists the protein pairs to run AF3 on: one row per unordered UniProt pair
(`uniprot_id_a` < `uniprot_id_b`) with at least one membership instance pair that has no structure.

- `reason`: comma-joined, from `no_af_metadata_entry`, `file_missing`, `file_unreadable` and
  `slice_failed`.
- Also `n_instance_pairs`, `labels` and `splits`.
- `all_methods.tsv` adds a `methods` column.

Nothing is hidden to equalise coverage. The task log only warns about the gaps.

### Reports

<details markdown="1">
<summary>Output files</summary>

- `reports/source_conflicts.tsv` — domain pairs dropped because two external sources disagreed on
  whether they interact.
- `reports/pruned_ddis.tsv` — DDIs dropped because a family has no domain instance to represent it.
- `reports/external_test_dropped.tsv` — external DDIs with no instance pair whose parent proteins
  differ.
- `reports/unresolved_split_rows.tsv` — split rows whose DDI or instance could not be resolved. Rows
  here mean a wiring error, not a data property.
- `reports/run_manifest.json` — what the run was made from: every param, the git commits of this
  repo and the ppi-splitting submodule, SHA-256 of the three manual inputs, `Last-Modified`/`ETag`
  of every download URL, and the positive/negative count of every (method, split). ILP solver
  status is not recorded: ppi-splitting does not publish it.
- `reports/split_qc/` (`REPORT_SPLIT_QC`, on the final master; every count is over
  `ddi_split_membership` rows, family/protein taken from the instance id, label from
  `domain_domain_interaction.negative`):
  - `split_overlap.tsv` — per method and unordered split pair, for each `unit` (`family`,
    `family_pair_positive`, `family_pair_negative`, `protein`, `instance`): the size of each split's
    set (`n_a`, `n_b`) and how many they share (`n_shared`).
  - `family_leaks.tsv` — `method, pfam_id, splits`: families that appear in more than one of
    `train`, `validation` and the test partition (all `test_*` splits together: the 1:1 and 1:10
    sets share their positives) of a method split by the ILP (every `split_method: 'ilp'` row of
    `workflows/domainsplit.nf`; the external `test` split is not part of that partition). The ILP
    guarantees there are none, so **any row fails the task** — and a failed task publishes nothing,
    so the leaks are also printed (first 50) to the task's stderr, which Nextflow shows in the error
    report. Other methods are reported in `split_overlap.tsv` only.
  - `test_variant_overlap.tsv` — per method with internal tests: its 1:1 test (`split_1_1`) vs
    `test_uniform_1_10`, per label, at `ddi_row` (DDI + instance pair) and `family_pair`
    granularity. Expected: positives identical, ILP negatives near-disjoint.
  - `external_test_profile.tsv` — per external-test method, its `test` split per label and source
    token (`__total__` = the whole label): family pairs, instance pairs, families, families also
    in the other label, families and family pairs also in train, and the min/median/max of the
    families' **train-positive degree** (distinct positive family pairs a family has in train)
    plus the median of each pair's smaller degree. High-degree families in one label only is the
    hub shortcut a model can score on without learning the interaction.
  - Only with `esm3_structure` in `--embedding_models_domainbench`:
    `esm3_structure_coverage.tsv` (per method, split, label: instance pairs with both instances in
    `esm3_structure_domain_embeddings.h5` vs total), `esm3_structure_missing_instances.tsv`
    (`instance_id, uniprot_id, reason, n_membership_rows`) and
    `esm3_structure_missing_accessions.tsv` (`accession, reason, n_instances`). `reason` is the
    protein's `afdb_structures_report.tsv` status (`not_in_afdb`, `sequence_mismatch`,
    `no_uniprot_sequence`, `download_failed`, `parse_failed`), `below_min_coverage` when the
    protein is `ok` (the domain fell under the structure-coverage or length cutoff), or
    `not_in_afdb_report`. Missing coverage only warns.
- `external_ddis/*.tsv` — each external source normalized, as the record of what it contributed
  before the merge/drop rules ran.

</details>

Every one of these is written even when empty (header only), so nothing downstream may block on a
report being non-empty.

### Splitting inputs and diagnostics

<details markdown="1">
<summary>Output files</summary>

- `ppi_splitting/inputs/families.txt` — the union of families instances were resolved for.
- `ppi_splitting/inputs/split_ddis.csv` — the exact positive set every ppi-splitting dataset row was
  given.
- `candidate_network/candidate_network.csv` — Pfam pairs the Y2H/MS screen tested without finding an
  interaction; the pool the `ilp_candidates` negatives are drawn from.
- `candidate_network/gene_pfam_mapping.json` — the gene → UniProt → Pfam lookup used to build it.
  Pass it back as `--negative_ppi_gene_mapping` to keep a re-run off the UniProt REST API.
- `dropped_families.tsv` — families Pfam had nothing usable for, with a `reason` column.
- `bias_analysis.html` — ppi-splitting's own MultiQC report on the splits it produced.

</details>

### Frozen splits

`frozen_splits/` (every run unless `--export_frozen_splits false`) lets a later run reuse exactly
this run's instances and splits, so AF3 predictions queued against this run stay valid:

- `fetch/`: FETCH_DOMAIN_META's `sequences.fasta`, `instances.tsv`, `species.tsv`,
  `go_annotations.tsv` and `dropped_families.tsv` (the sampled domain instances).
- `split_ddis/split_ddis.csv`: the positive population handed to ppi-splitting.
- `splits/{family,instance}/<method>/<split>.csv`: what ppi-splitting handed to INGEST_SPLITS.
- `bundle_manifest.json`: SHA-256 of every file, the split list, the params the bundle depends on
  (`frozenParamKeys()` in `workflows/domainsplit.nf`), and the final master's membership and
  protein-sequence fingerprints.

`--frozen_splits <dir>` skips FETCH_DOMAIN_META and PPI_SPLITTING. Three checks guard it:

1. At launch: the params and the (method, split) set must match the bundle.
2. VERIFY_FROZEN_SPLITS, before ingest: every bundle file hashes to its entry, and this run's
   positive population (rebuilt from its own sources) equals the frozen one.
3. CHECK_FROZEN_DB, before structure and the split databases: the final master's membership rows
   (external test included) and protein sequences match the frozen run's fingerprints.

Pin every source (local copies for each `url_*` and `uniprot_dat_urls`) to the frozen run's
files, or check 2 or 3 fails. A frozen run does not republish ppi-splitting's own reports
(`union/data/dropped_families.tsv`, the bias report): those belong to the original run.
`reports/run_manifest.json` records both fingerprints (`membership_fingerprint`,
`protein_fingerprint`) for every run.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
