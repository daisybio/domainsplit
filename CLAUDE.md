# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

`daisybio/domainsplit` is an nf-core-flavored Nextflow DSL2 pipeline (template v4.0.2, but `is_nfcore: false`) that builds a SQLite database of domain-domain interactions (`domainsplit.sqlite3`) from public + manually-supplied sources (3did, HIPPIE, PPIDM, Negatome, a Y2H/MS negative-PPI screen, UniProt, STRING, Pfam, pfam2go), enriches it with ProtT5 + ESM3/ESMC embeddings and GO terms, then carves it into one database per (splitting strategy, split).

Splitting and negative sampling are **delegated to `ppi-splitting-pipeline`**, imported as a git submodule at `subworkflows/external/ppi-splitting` and called as the `PPI_SPLITTING` named workflow with three dataset rows (`random`, `minimal_leakage`, `external_test`). `PLAN_ppi_splitting_integration.md` is the spec for that integration and the record of what is and is not wired; read it before changing anything in that path.

## Running the pipeline

Standard invocations (`nextflow run . -profile test,docker`, `nf-test test`, `nf-core pipelines lint`, `pre-commit run --all-files`, `nf-core modules install`) work as expected. Two non-obvious rules:

- **Pass pipeline parameters via CLI flags or `-params-file` only.** Per nf-core convention, parameters MUST NOT be set in `-c` config files — only process/executor config goes there.
- `conf/test.config` fixtures are meant for `-stub` runs (see `tests/default.nf.test`): process bodies never execute, so those tests verify DAG wiring (channel topology, split fan-out, publish paths) and nothing about per-module logic. **That test currently cannot run**: `-stub` silently falls back to the real `script:` block for a process with no `stub:` block, and no ppi-splitting process defines one. `nextflow run . -profile test -preview` is the working DAG check meanwhile.
- Add `NXF_OFFLINE=1` locally, or the `nfcore_custom.config` fetch times out and config parsing fails.
- `conf/test_full.config` is still an unmodified nf-core template placeholder — not usable as-is.

## Gotchas

- **Everything is a value channel** (a single SQLite path) flowing between processes, NOT a per-sample samplesheet. The `samplesheet` channel from `PIPELINE_INITIALISATION` is wired in but `DOMAINSPLIT` ignores it — input comes entirely from `params.url_*`. Don't reshape `domainsplit_db` channels into meta-tuple channels without updating every downstream consumer.
- Three required inputs have **no download automation and no documented provenance/versioning**: `hippie_tsv`, `ppidm_tsv`, `negative_ppi_parquet`. They must be supplied manually; no `url_*` param exists for them.
- Split DBs publish to `databases/<method>/<split>.sqlite3` (see the `output:` block in `main.nf`). The README claims `split_databases/...` and is out of date. The README also still carries unfilled nf-core template TODOs (intro summary, workflow figure, step bullets).
- `HF_TOKEN` (HuggingFace, for ESM3/ESMC weights) is required via `.env` / `nextflow secrets` — see `.env.example`.
- Some processes shell out to `bin/mysql2sqlite` (vendored awk script) — keep `bin/` executable.
- Test coverage is thin: a handful of pytest files in `tests/python/` over individual `bin/` scripts, plus stub-only nf-test DAG wiring. Most modules have no unit tests.

## Conventions

- Module imports: `include { NAME } from '../modules/local/<dir>/main.nf'`. Local subworkflow files are referenced via `./subworkflows/...` **relative to `workflows/domainsplit.nf`, NOT the project root**.
- nf-core modules are installed via `nf-core modules install` and tracked in `modules.json`; don't hand-edit them — override behavior with `ext.args` in `conf/modules.config`. Only nf-core _subworkflows_ are installed at the moment; `modules/nf-core/` is empty.
- Never edit anything under `subworkflows/external/ppi-splitting/` — it is a submodule of another repository. Nextflow reads only the root project's `nextflow.config`, so that pipeline's params come in via `includeConfig` of its `conf/params.config` and everything else it needs (container, per-process resources, the `gurobi` label, publishing) is restated in `conf/ppi_splitting.config`. Its `bin/` is baked into its Docker image, so the image tag and the submodule pin move together.
- `conf/modules.config` is also where `publishDir` lives; the default path tokenizes the process name (`${process.split(':').last().split('_').first().lowercase()}`) into `${params.outdir}/<group>`.
- Run `pre-commit install` after cloning (`.pre-commit-config.yaml` runs prettier + editorconfig-checker + nf-core linting).
- `.nf-core.yml` lists files intentionally diverging from the nf-core template (under `lint.files_unchanged` / `lint.files_exist`) — adding files in those paths re-enables the lint check.
- `is_nfcore: false` — this pipeline borrows the template but is not under the `nf-core/` org, so don't add `nf-core/` paths or AWS test workflows.

<!-- code-review-graph MCP tools -->

## MCP Tools: code-review-graph

**IMPORTANT: This project has a knowledge graph. ALWAYS use the code-review-graph MCP
tools BEFORE using Grep/Glob/Read to explore the codebase.** The graph is faster, cheaper
(fewer tokens), and gives structural context (callers, dependents, test coverage) that
file scanning cannot. The graph auto-updates on file changes via hooks.

- **Exploring code**: `semantic_search_nodes` or `query_graph` instead of Grep
- **Understanding impact**: `get_impact_radius` / `get_affected_flows` instead of manually tracing imports
- **Code review**: `detect_changes` + `get_review_context` instead of reading entire files
- **Finding relationships**: `query_graph` with `callers_of` / `callees_of` / `imports_of` / `tests_for`
- **Architecture questions**: `get_architecture_overview` + `list_communities`
- **Renames / dead code**: `refactor_tool`

Fall back to Grep/Glob/Read **only** when the graph doesn't cover what you need.
