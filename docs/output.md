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

</details>

With the default `--ppi_splitting_multi_negset true` there are **5 method directories and 18
databases**:

| method                 | splits                                                   |
| ---------------------- | -------------------------------------------------------- |
| `random`               | `train`, `validation`, `test_balanced`, `test_realistic` |
| `minimal_leakage`      | `train`, `validation`, `test_balanced`, `test_realistic` |
| `minimal_leakage_hcni` | `train`, `validation`, `test_balanced`, `test_realistic` |
| `external_test`        | `train`, `validation`, `test`                            |
| `external_test_hcni`   | `train`, `validation`, `test`                            |

`random` is the naive baseline (random split, uniform negatives). `minimal_leakage` is the main
product (ILP split, ILP negatives). The `_hcni` directories are the same positive split with
**h**igh-**c**onfidence **n**egative **i**nstances: negatives drawn from the candidate network of
Pfam pairs a Y2H/MS screen tested without finding an interaction. `external_test`'s `test` split is
not an internal held-out slice at all — it is the external DDIs (PPIDM / single-domain-PPI /
Negatome), unseen by construction.

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
h5[str(domain_id)][instance_id]
```

where `domain_id` is `domain.id` and `instance_id` is `domain_protein_map.instance_id` — precisely
`h5[str(domain_id)][COALESCE(instance_id, 'r' || rowid)]` over the split database's
`domain_protein_map`.

Only the **cut domain sequence** is embedded, never the parent protein. A per-residue protein
embedding sliced to a domain's coordinates still carries protein context, and protein context
correlates with the interaction partner — which would smuggle protein identity into a split that
was partitioned on families.

Root attributes on each file: `model`, `pooling`, `dim`, `dtype`, `key_layout`, `n_domains`,
`n_instances`, `domainsplit_run`.

> [!IMPORTANT]
> `domain.id` is a **surrogate integer**. One embedding file is valid across every split database of
> the same run, and silently wrong across runs. Check `domainsplit_run` before pairing an embedding
> file with databases you did not produce together.

Which models run is `--embedding_models_domainbench` (default `esm3,esmc,prott5`); set it empty to
skip embedding entirely. ESM3 and ESMC need `HF_TOKEN`; ProtT5 does not.

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

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
