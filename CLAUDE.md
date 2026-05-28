# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

`daisybio/domainsplit` is an nf-core-flavored Nextflow DSL2 pipeline (template v4.0.2, but `is_nfcore: false`) that builds a SQLite database of domain-domain interactions (`domainsplit.sqlite3`) from public sources (3did, UniProt, Negatome, STRING, Pfam, pfam2go) and then splits it into train/validation/test partitions using several strategies, with leakage reduction as the main scientific goal.

## Common commands

Run pipeline (test profile):

```bash
nextflow run . -profile test,docker --outdir results
```

Run pipeline (real run, override URLs/paths via params file):

```bash
nextflow run . -profile docker --outdir results -params-file params.yml
```

Run a single nf-test test file:

```bash
nf-test test tests/default.nf.test
```

Run all nf-tests (uses `profile = "test"` from `nf-test.config`):

```bash
nf-test test
```

Lint (matches the `linting.yml` CI workflow):

```bash
nf-core pipelines lint
pre-commit run --all-files
prettier --check .
```

Install / update an nf-core module:

```bash
nf-core modules install <tool/subtool>
nf-core modules update <tool/subtool>
```

> Pass pipeline parameters via CLI flags or `-params-file`. Per nf-core convention, parameters MUST NOT be set in `-c` config files — only process/executor config goes there.

## Architecture

Entrypoint chain:

- `main.nf` — runs `PIPELINE_INITIALISATION` → `DAISYBIO_DOMAINSPLIT` (wraps `DOMAINSPLIT`) → `PIPELINE_COMPLETION`.
- `workflows/domainsplit.nf` — top-level scientific workflow. Reads all input URLs/files from `params.url_*`, chains the subworkflows below, and declares a `publish:`/`output:` block that routes split databases into `split_databases/${method}/${split}.sqlite3`.
- `modules/local/init_domainsplit_db/` (`INIT_DOMAINSPLIT_DB`) — creates the empty `domainsplit.sqlite3` schema that flows through every later stage.
- `subworkflows/local/collect_ddi_data/` (`COLLECT_DDI_DATA`) — downloads 3did + Negatome, inserts positive/negative DDIs, applies optional smoke filter. `sqlite_3did` is internal and does NOT escape this subworkflow.
- `subworkflows/local/curate_domains/` (`CURATE_DOMAINS`) — extracts unique Pfam IDs from the in-build DB's `domain_domain_interaction` table (inline sqlite3, no python), downloads Pfam alignments, and creates the protein↔domain map.
- `subworkflows/local/generate_embeddings/` (`GENERATE_EMBEDDINGS`) — parallel ProtT5 (per-residue HDF5) + ESM3/ESMC (per-residue protein + pooled domain) embedding generation.
- `subworkflows/local/enrich_ddi_database/` (`ENRICH_DDI_DATABASE`) — sequential chain of five `INSERT_*` processes (`INSERT_DOMAIN_GO_TERMS`, `INSERT_PROTEINS_WITH_EMBEDDINGS`, `INSERT_PROTEIN_GO_TERMS`, `INSERT_PPI`, `INSERT_DOMAIN_PROTEIN_MAPPING`). Each opens the SQLite emitted by the previous step, performs one phase, commits, and emits the DB forward.
- `subworkflows/local/split_domainsplit_database/` (`SPLIT_DOMAINSPLIT_DATABASE`) — extracts domain sequences, clusters with `MMSEQS_EASYCLUSTER` (nf-core module), and runs split strategies producing per-split SQLite DBs: `RANDOM_DDI_SPLIT` (biased baseline) and `MINIMAL_LEAKAGE_SPLIT_DOMAIN` (spectral graph-partitioning). `map_split_files()` is the helper that re-keys flattened split outputs into `[meta, path]` tuples where `meta = [id, split, method]`.
- `subworkflows/local/utils_nfcore_domainsplit_pipeline/` — pipeline init/completion/methods-description helpers (template-generated).

Modules:

- `modules/local/*` — all the scientific work: 3did SQL→SQLite via `bin/mysql2sqlite`, DDI insert/smoke-filter, Pfam alignment + protein-domain mapping, the five `enrich/insert_*` phase modules (each calling a corresponding `bin/insert_*.py` script), ProtT5/ESM embedding generation, the splitter modules, sequence extractors.
- `modules/nf-core/mmseqs/easycluster/` — only nf-core module currently installed. Pinned in `modules.json`.

Config layout:

- `nextflow.config` — global `params` block (all the `url_*` defaults live here), profile definitions (`docker`, `singularity`, `conda`, `podman`, `apptainer`, `charliecloud`, `test`, `test_full`, `arm`, `debug`), and registry settings (`quay.io` for all container engines).
- `conf/modules.config` — per-process `publishDir` and `ext.args`. Default publish path tokenizes the process name (`${process.split(':').last().split('_').first().lowercase()}`) into `${params.outdir}/<group>`. Tweak `ext.args` here, NOT inside module files — example: `MMSEQS_EASYCLUSTER` clustering thresholds (`--min-seq-id 0.4 -c 0.8`).
- `conf/base.config` — resource labels (`process_low/medium/high/long`).
- `conf/test.config`, `conf/test_full.config` — test profiles. The current `input` value is still the nf-core template placeholder (`viralrecon` samplesheet); update before relying on the `test` profile.
- `nextflow_schema.json` — parameter schema; validated when `params.validate_params = true`.

Data flow shape: everything is a value channel (single SQLite path) flowing between processes, NOT a per-sample samplesheet. The `samplesheet` channel from `PIPELINE_INITIALISATION` is wired in but `DOMAINSPLIT` ignores it — input comes entirely from `params.url_*`. Keep this in mind when editing: don't reshape `domainsplit_db` channels into meta-tuple channels without updating every downstream consumer.

## Conventions

- DSL2 only. Module imports use `include { NAME } from '../modules/local/<dir>/main.nf'`. Local subworkflow files live alongside `workflows/domainsplit.nf` and are referenced via `./subworkflows/...` from within that file (relative to the workflow file, NOT the project root).
- nf-core modules are installed via `nf-core modules install` and tracked in `modules.json`; don't hand-edit them — override behavior with `ext.args` in `conf/modules.config`.
- Pre-commit hooks (`.pre-commit-config.yaml`) run prettier + editorconfig-checker + nf-core linting. Run `pre-commit install` after cloning.
- `.nf-core.yml` lists files intentionally diverging from the nf-core template (under `lint.files_unchanged` / `lint.files_exist`) — adding files in those paths will re-enable the lint check.
- `is_nfcore: false` in `.nf-core.yml` — this pipeline borrows the template but is not under the `nf-core/` org, so don't add `nf-core/` paths or AWS test workflows.

## Notes from README

- Main output is `domainsplit.sqlite3` plus per-method split DBs in `split_databases/<method>/<split>.sqlite3`.
- Some processes shell out to `bin/mysql2sqlite` (vendored awk script) — keep `bin/` executable.

<!-- code-review-graph MCP tools -->

## MCP Tools: code-review-graph

**IMPORTANT: This project has a knowledge graph. ALWAYS use the
code-review-graph MCP tools BEFORE using Grep/Glob/Read to explore
the codebase.** The graph is faster, cheaper (fewer tokens), and gives
you structural context (callers, dependents, test coverage) that file
scanning cannot.

### When to use graph tools FIRST

- **Exploring code**: `semantic_search_nodes` or `query_graph` instead of Grep
- **Understanding impact**: `get_impact_radius` instead of manually tracing imports
- **Code review**: `detect_changes` + `get_review_context` instead of reading entire files
- **Finding relationships**: `query_graph` with callers_of/callees_of/imports_of/tests_for
- **Architecture questions**: `get_architecture_overview` + `list_communities`

Fall back to Grep/Glob/Read **only** when the graph doesn't cover what you need.

### Key Tools

| Tool                        | Use when                                               |
| --------------------------- | ------------------------------------------------------ |
| `detect_changes`            | Reviewing code changes — gives risk-scored analysis    |
| `get_review_context`        | Need source snippets for review — token-efficient      |
| `get_impact_radius`         | Understanding blast radius of a change                 |
| `get_affected_flows`        | Finding which execution paths are impacted             |
| `query_graph`               | Tracing callers, callees, imports, tests, dependencies |
| `semantic_search_nodes`     | Finding functions/classes by name or keyword           |
| `get_architecture_overview` | Understanding high-level codebase structure            |
| `refactor_tool`             | Planning renames, finding dead code                    |

### Workflow

1. The graph auto-updates on file changes (via hooks).
2. Use `detect_changes` for code review.
3. Use `get_affected_flows` to understand impact.
4. Use `query_graph` pattern="tests_for" to check coverage.
