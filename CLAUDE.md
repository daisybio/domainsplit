# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

`daisybio/domainsplit` is an nf-core-flavored Nextflow DSL2 pipeline (template v4.0.2, but `is_nfcore: false`) that builds a SQLite database of domain-domain interactions (`domainsplit.sqlite3`) from public + manually-supplied sources (3did, HIPPIE, PPIDM, Negatome, a Y2H/MS negative-PPI screen, UniProt, STRING, Pfam, pfam2go), enriches it with ProtT5 + ESM3/ESMC embeddings and GO terms, then splits it into train/validation/(test|optimization) partitions using several leakage-reduction strategies, run once per negative-DDI sampling method (`deletion` / `random_addition`).

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

- `main.nf` — runs `PIPELINE_INITIALISATION` → `DAISYBIO_DOMAINSPLIT` (wraps `DOMAINSPLIT`) → `PIPELINE_COMPLETION`. Its `publish:`/`output:` block routes `domainsplit_db` and `split_db` (→ `databases/${method}/${split}.sqlite3` — NOT `split_databases/`, despite what the README says). `DOMAINSPLIT`'s `bias_report` output is **not** wired into this block or into `DAISYBIO_DOMAINSPLIT`'s emit — `ANALYZE_DDI_BIAS` runs but its report is never published anywhere.
- `workflows/domainsplit.nf` — top-level scientific workflow. Reads all input URLs/files from `params.url_*`, plus three required non-URL file params with no download automation (`hippie_tsv`, `ppidm_tsv`, `negative_ppi_parquet` — must be supplied manually), chains the subworkflows below.
- `modules/local/init_domainsplit_db/` (`INIT_DOMAINSPLIT_DB`) — creates the empty `domainsplit.sqlite3` schema that flows through every later stage.
- `subworkflows/local/collect_ddi_data/` (`COLLECT_DDI_DATA`) — inserts every DDI source in a fixed order, so the earlier source wins on a duplicate domain pair (`INSERT OR IGNORE`):
  1. `INSERT_3DID` — 3did positives (raw SQLite via `DOWNLOAD_3DID_SQLITE` + `bin/mysql2sqlite`). `sqlite_3did` is internal and does NOT escape this subworkflow.
  2. `INSERT_SINGLE_DOMAIN_PPI` — positives inferred from HIPPIE PPIs between two single-domain proteins, using a reviewed-human SwissProt→Pfam map (`BUILD_SWISSPROT_PFAM_MAP`).
  3. `INSERT_PPIDM` — PPIDM predicted positives, source-tagged by class (`PPIDM_Bronze/Silver/Gold`).
  4. `INSERT_NEGATOME` — Negatome negatives.
  5. optional `REMOVE_SELF_INTERACTIONS` — only runs when `params.self_interaction = false` (default `true`, so this is skipped by default).
  6. High-confidence non-PPI negatives from a Y2H/MS screen parquet, via uncapped degree-proportional node-pair sampling (DANS): `BUILD_PPI_NEGATIVE_POOL` maps bait/prey genes → UniProt → Pfam (REST) once, then `SELECT_PPI_NEGATIVE_DANS` runs **twice** — once per negative-DDI method, `deletion` and `random_addition` — and `INSERT_PPI_NEGATIVE_SELECTION` writes both under distinct source labels (`3did_<method>` / `inferred_ppi_screen_negative_for_<method>`).
  7. optional `SMOKE_FILTER` if `params.smoke_test_n_ddis` is set.
- `subworkflows/local/curate_domains/` (`CURATE_DOMAINS`) — extracts unique Pfam IDs from the in-build DB's `domain_domain_interaction` table (inline sqlite3, no python), downloads Pfam alignments in batches (`params.pfam_download_batch_size`), and creates the protein↔domain map.
- `modules/local/esm_embeddings/` (`generate_esm_embeddings`) — called directly from `workflows/domainsplit.nf` (not a subworkflow). ESM3/ESMC per-residue protein (sharded) + pooled domain embeddings, generated in parallel with the ProtT5 path consumed later in `ENRICH_DDI_DATABASE`. Requires `HF_TOKEN` (see `.env.example`).
- `subworkflows/local/enrich_ddi_database/` (`ENRICH_DDI_DATABASE`) — sequential chain of five `INSERT_*` processes (`INSERT_DOMAIN_GO_TERMS`, `INSERT_PROTEINS_WITH_EMBEDDINGS` [ProtT5 HDF5 + ESM protein embeddings], `INSERT_PROTEIN_GO_TERMS`, `INSERT_PPI` [STRING], `INSERT_DOMAIN_PROTEIN_MAPPING` [+ ESM domain embeddings]). Each opens the SQLite emitted by the previous step, performs one phase, commits, and emits the DB forward.
- `modules/local/analyze_ddi_bias/` (`ANALYZE_DDI_BIAS`) — runs once on the fully-enriched DB, emits a `bias_analysis/` report dir. See the publish-gap note on `main.nf` above.
- `subworkflows/local/split_domainsplit_database/` (`SPLIT_DOMAINSPLIT_DATABASE`) — extracts domain sequences and clusters them once with `MMSEQS_EASYCLUSTER` (nf-core module); every leakage-aware partition below reuses the same clusters. Each strategy runs **once per negative-DDI method** (`deletion` / `random_addition`), yielding 6 output folders:
  - `random_ddi_<method>` — `RANDOM_DDI_SPLIT`, biased baseline (random 60/20/20 train/optimization/test).
  - `minimal_leakage_domain_<method>` — `MINIMAL_LEAKAGE_SPLIT_DOMAIN`, spectral graph-partitioning (Laplacian eigenvectors + weighted k-means + Kernighan-Lin refinement) over MMseqs2 domain clusters, same 60/20/20 fractions, over all DDI sources.
  - `external_validation_<method>` — `MINIMAL_LEAKAGE_SPLIT_DOMAIN` again, but with its `source_filter` restricted to that method's "core" sources only (its 3did copy + its own PPI-screen negatives), 80/20 train/validation. The `test` split is built separately by `SUBSET_DDIS_BY_SOURCE`, which keeps the complementary held-out sources (`single_domain_ppi`, `PPIDM_Bronze/Silver/Gold`, `negatome`) as-is (no leakage-aware partitioning, no per-method distinction) — the same test DB is routed into both `external_validation_deletion/test` and `external_validation_random_addition/test`. Net effect: train/validate with leakage control on one curated-source population, then test generalization against DDIs from entirely different, independently-curated source populations.
  - `map_split_dbs()` is the helper that re-keys flattened split outputs into `[meta, path]` tuples where `meta = [id: "${method}_${split}", split, method]`.
- `subworkflows/local/utils_nfcore_domainsplit_pipeline/` — pipeline init/completion/methods-description helpers (template-generated).

Modules:

- `modules/local/*` — all the scientific work: 3did SQL→SQLite via `bin/mysql2sqlite`, one download/insert module pair per DDI source (see `collect_ddi_data` above), Pfam alignment + protein-domain mapping, the five `enrich/insert_*` phase modules (each calling a corresponding `bin/insert_*.py` script), ProtT5/ESM embedding generation, the bias-analysis module, the splitter modules, sequence extractors.
- `modules/nf-core/mmseqs/easycluster/` — only nf-core module currently installed. Pinned in `modules.json`.

Config layout:

- `nextflow.config` — global `params` block (all the `url_*` defaults live here), profile definitions (`docker`, `singularity`, `conda`, `podman`, `apptainer`, `charliecloud`, `test`, `test_full`, `arm`, `debug`), and registry settings (`quay.io` for all container engines).
- `conf/modules.config` — per-process `publishDir` and `ext.args`. Default publish path tokenizes the process name (`${process.split(':').last().split('_').first().lowercase()}`) into `${params.outdir}/<group>`. Tweak `ext.args` here, NOT inside module files — example: `MMSEQS_EASYCLUSTER` clustering thresholds (`--min-seq-id 0.4 -c 0.8`).
- `conf/base.config` — resource labels (`process_low/medium/high/long`).
- `conf/test.config` — points every `url_*`/file param at tiny local fixtures under `tests/data/`; run with `-stub` (see `tests/default.nf.test`), so process content is never read, only DAG wiring (channel topology, split fan-out, publish paths) is verified.
- `conf/test_full.config` — still an unmodified nf-core template placeholder (`input`/`fasta` = `pipelines_testdata_base_path`); not usable as-is.
- `nextflow_schema.json` — parameter schema; validated when `params.validate_params = true`.

Data flow shape: everything is a value channel (single SQLite path) flowing between processes, NOT a per-sample samplesheet. The `samplesheet` channel from `PIPELINE_INITIALISATION` is wired in but `DOMAINSPLIT` ignores it — input comes entirely from `params.url_*`. Keep this in mind when editing: don't reshape `domainsplit_db` channels into meta-tuple channels without updating every downstream consumer.

## Conventions

- DSL2 only. Module imports use `include { NAME } from '../modules/local/<dir>/main.nf'`. Local subworkflow files live alongside `workflows/domainsplit.nf` and are referenced via `./subworkflows/...` from within that file (relative to the workflow file, NOT the project root).
- nf-core modules are installed via `nf-core modules install` and tracked in `modules.json`; don't hand-edit them — override behavior with `ext.args` in `conf/modules.config`.
- Pre-commit hooks (`.pre-commit-config.yaml`) run prettier + editorconfig-checker + nf-core linting. Run `pre-commit install` after cloning.
- `.nf-core.yml` lists files intentionally diverging from the nf-core template (under `lint.files_unchanged` / `lint.files_exist`) — adding files in those paths will re-enable the lint check.
- `is_nfcore: false` in `.nf-core.yml` — this pipeline borrows the template but is not under the `nf-core/` org, so don't add `nf-core/` paths or AWS test workflows.

## Notes from README

- Main output is `domainsplit.sqlite3` plus per-method split DBs. README says `split_databases/<method>/<split>.sqlite3`; actual publish path (`main.nf` output block) is `databases/<method>/<split>.sqlite3` — README is out of date on this.
- README still has unfilled nf-core template TODOs (intro summary, workflow figure, step bullets) and never documents `hippie_tsv`/`ppidm_tsv`/`negative_ppi_parquet` — no URL params exist for these, they're pure manual acquisition with no documented provenance/versioning.
- Some processes shell out to `bin/mysql2sqlite` (vendored awk script) — keep `bin/` executable.
- `HF_TOKEN` (HuggingFace, for ESM3/ESMC weights) is required via `.env`/`nextflow secrets` — see `.env.example`.
- Test coverage is thin: `tests/python/` only covers `insert_negatome`, `insert_ppidm`, `insert_ppi_negative_selection`, `select_ppi_negative_dans`. No unit tests for `insert_3did`, `build_ppi_negative_pool`, `build_swissprot_pfam_map`, the GO/PPI/embedding inserts, `analyze_ddi_bias`, or any split module — nf-test only exercises stub-run wiring, not per-module logic.

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
