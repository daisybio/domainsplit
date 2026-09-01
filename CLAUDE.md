# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

`daisybio/domainsplit` is an nf-core-flavored Nextflow DSL2 pipeline (template v4.0.2, but `is_nfcore: false`) that builds a SQLite database of domain-domain interactions (`domainsplit.sqlite3`) from public + manually-supplied sources (3did, HIPPIE, PPIDM, Negatome, a Y2H/MS negative-PPI screen, UniProt, STRING, Pfam, pfam2go), enriches it with GO terms, protein sequences and STRING PPIs, then carves it into one database per (splitting strategy, split). Embeddings are **not** in the database: every domain instance is embedded with ProtT5 + ESM3/ESMC and published as three HDF5 files under `embeddings/`. The master carries no BLOBs.

Splitting and negative sampling are **delegated to `ppi-splitting-pipeline`**, imported as a git submodule at `subworkflows/external/ppi-splitting` and called as the `PPI_SPLITTING` named workflow with three dataset rows (`random`, `minimal_leakage`, `external_test`). The two ILP rows each ask for two negative sets (`ilp` + `ilp_candidates`), so three rows yield **5 method directories / 18 split databases**. `PLAN_ppi_splitting_integration.md` is the spec for that integration and the record of what is and is not wired; read it before changing anything in that path.

## Running the pipeline

Standard invocations (`nextflow run . -profile test,docker`, `nf-test test`, `nf-core pipelines lint`, `pre-commit run --all-files`, `nf-core modules install`) work as expected. Two non-obvious rules:

- **Pass pipeline parameters via CLI flags or `-params-file` only.** Per nf-core convention, parameters MUST NOT be set in `-c` config files — only process/executor config goes there.
- `conf/test.config` fixtures are meant for `-stub` runs (see `tests/default.nf.test`): process bodies never execute, so those tests verify DAG wiring (channel topology, split fan-out, publish paths) and nothing about per-module logic. **That test currently cannot run**: `-stub` silently falls back to the real `script:` block for a process with no `stub:` block, and no ppi-splitting process defines one. `nextflow run . -profile test -preview` is the working DAG check meanwhile.
- Add `NXF_OFFLINE=1` locally, or the `nfcore_custom.config` fetch times out and config parsing fails.
- `conf/test_full.config` is still an unmodified nf-core template placeholder — not usable as-is.

## Gotchas

- **Everything is a value channel** (a single SQLite path) flowing between processes, NOT a per-sample samplesheet. The `samplesheet` channel from `PIPELINE_INITIALISATION` is wired in but `DOMAINSPLIT` ignores it — input comes entirely from `params.url_*`. Don't reshape `domainsplit_db` channels into meta-tuple channels without updating every downstream consumer.
- Three required inputs have **no download automation and no documented provenance/versioning**: `hippie_tsv`, `ppidm_tsv`, `negative_ppi_parquet`. They must be supplied manually; no `url_*` param exists for them.
- Split DBs publish to `databases/<method>/<split>.sqlite3` (see the `output:` block in `main.nf`); domain embeddings publish to `embeddings/<model>_domain_embeddings.h5` via a `withName` rule in `conf/modules.config`. An `EXPORT_*`/`GENERATE_*` process without such a rule publishes to a directory tokenized from its own name, which is how a mis-routed file appears. The README still carries unfilled nf-core template TODOs (intro summary, workflow figure, step bullets).
- `HF_TOKEN` (HuggingFace, for the gated ESM3/ESMC weights) is required via `.env` / `nextflow secrets` — see `.env.example`. ProtT5 is ungated and needs no token.
- The UniProt flat files in `uniprot_dat_urls` feed **five** consumers: GO terms, the accession→Pfam map, the protein FASTA, the STRING-id map (`DR STRING;`, which replaced the per-organism `url_uniprot_id_mapping` download) and `FETCH_DOMAIN_META`'s protein universe. `--dat` is repeatable and **Swiss-Prot must be listed first** — an accession promoted from TrEMBL to Swiss-Prot between releases is in both files and both parsers take the first writer. The accession→Pfam map is written for **reviewed entries only** (`--pfam-map-reviewed-only`, default on): it feeds `parse_single_domain_ppi.py`, which infers a DDI from a HIPPIE PPI between two _single-domain_ proteins, and on a TrEMBL entry "exactly one domain" is not a claim an automatic annotation over a possible fragment can support. HIPPIE is read by **UniProt accession only**; a file without `uniprot_accession_A/B` columns is rejected rather than silently yielding zero DDIs.
- `domain.id` is a surrogate integer, so **nothing published keys on it**. The embedding HDF5s key on the Pfam accession (`h5[pfam_id][instance_id]`, `domain` is `UNIQUE(pfam_id)` so the two are 1:1), and carry no run identifier — accession keys need no cross-run guard. Anything new that leaves this pipeline must key on `pfam_id` too — a surrogate id crossing the boundary resolves to the wrong row in another run instead of failing.
- **`instance_tier` is the protein universe, and three params are derived from it.** One knob, four alternatives (not a ladder): `human_reviewed` (default, ~20 k entries), `all_species_reviewed` (~573 k), `human_any_review_status` (~200 k) and `all_species_any_review_status`, which **hard-fails** — it needs the 110 GB TrEMBL file and `FETCH_DOMAIN_META` parses the universe into an in-memory dict (daisybio/domainsplit#4). From it come `uniprot_dat_urls` (which files are downloaded), `swissprot_taxon_ids` (which taxa survive parsing) and `instance_tiers` (ppi-splitting's stratum-name list).
  - **The derivation lives in `nextflow.config`, not in the workflow**, because `FETCH_DOMAIN_META` is in the read-only submodule and reads `params.instance_tiers` from its own script block, and `params` is read-only by the time a workflow body runs. A CLI `--instance_tier X` _is_ visible to the config — Nextflow merges command-line params into the binding before parsing it — which is what makes that work; each derived value guards on `!= null`, not truthiness, so an explicitly empty value survives. `tests/python/test_instance_tier_wiring.py` pins this.
  - `instance_tier` → `instance_tiers` is **the only place the two vocabularies meet**. Upstream has no notion of `human_any_review_status`; this repo has no notion of `other_reviewed`.
  - Upstream fills strata in the order `human_reviewed → other_reviewed → human_unreviewed → other_unreviewed`: **reviewed outranks human**. Lower strata fill _freely_, not as a top-up, so widening the universe reshapes the instance pool of families that were never dropped — only `human_reviewed` reproduces the previous run's numbers.
  - The retired vocabulary is `human_only` (→ `human_reviewed`) and `any` (→ `all_species_reviewed`, **not** `all_species_any_review_status`). Both are rejected with the mapping rather than aliased, because aliasing `any` would silently change the universe of a working invocation.
  - `ingest_instances.py` asserts the taxon invariant at the other end; `protein.reviewed` (from `instances.tsv`'s `source_db`) and `protein.taxon_id` make the result stratifiable, and `reports/ddi_tier_breakdown.tsv` counts surviving DDIs per stratum — a DDI above `human_reviewed` is one a narrower universe would have pruned.
- **Domain instances come from `Pfam-A.regions.tsv.gz`, never `Pfam-A.fasta`.** Pfam publishes the FASTA 90 % non-redundant, which drops the human member of every well-conserved family — a 12,769-family human request kept only 2,140. The regions file has no redundancy reduction; it carries coordinates only, so the UniProt flat files supply sequence (`SQ`), taxon (`OX`) and the review flag (the `ID` line's `Reviewed;`/`Unreviewed;` token). That fix lives in the submodule's `bin/fetch_domains.py`.
- **STRING is fetched per organism, never as the monolith.** `protein.links.v12.0.txt.gz` (all organisms) is 138 GB and >99 % irrelevant, so `FETCH_STRING_LINKS` downloads `<taxid>.protein.links.v<rel>.txt.gz` per taxon and drops every edge whose endpoints are not both in this run's `protein` table — `insert_ppi.py` therefore sees a small file and streams it. **The taxon list comes from the STRING id prefixes (`9606.ENSP…`), not from `protein.taxon_id`**: UniProt's `OX` can be strain-level with no STRING file, while the prefix is STRING's own species id, so the join is exact by construction. A 404 per taxon is counted in `reports/string_taxa_report.tsv`, never fatal; _every_ taxon 404ing is. Setting `url_string` pins one explicit links file and skips the fetch, which is how `-profile test` stays offline. **TrEMBL proteins get no PPI enrichment and no links file changes that** — UniProt does not cross-reference STRING for unreviewed entries (measured: 4,616 `uniprot_trembl_human` entries carried 15 `DR STRING;` lines), so both parsers report the gap rather than leaving it to be re-investigated.
- **One cache root: `cache_dir`.** ppi-splitting's Pfam/UniProt downloads land under `pfam-<release>/` and the per-organism STRING files under `string/`. `interpro_cache` is ppi-splitting's own name for it and is kept as a legacy alias — `nextflow.config` sets one from the other, `cache_dir` wins.
- **A boolean param cannot be tested for truthiness directly.** A CLI `--flag false` arrives as the String `"false"`, which is truthy in Groovy, so `params.flag ? a : b` silently takes the `a` branch. Use `params.flag.toString().toBoolean()`. This bit `ppi_splitting_multi_negset` and `embedding_require_gpu`, both of which are documented off-switches.
- Some processes shell out to `bin/mysql2sqlite` (vendored awk script) — keep `bin/` executable.
- Test coverage is thin: a handful of pytest files in `tests/python/` over individual `bin/` scripts, plus stub-only nf-test DAG wiring. Most modules have no unit tests. `tests/python/test_instance_tier_wiring.py` is the exception that shells out to `nextflow run -preview`; it skips when `nextflow` is not on PATH.

## Conventions

- Module imports: `include { NAME } from '../modules/local/<dir>/main.nf'`. Local subworkflow files are referenced via `./subworkflows/...` **relative to `workflows/domainsplit.nf`, NOT the project root**.
- nf-core modules are installed via `nf-core modules install` and tracked in `modules.json`; don't hand-edit them — override behavior with `ext.args` in `conf/modules.config`. Only nf-core _subworkflows_ are installed at the moment; `modules/nf-core/` is empty.
- **`subworkflows/external/ppi-splitting/` is read-only. The only permitted operation is pulling the current pushed version** (`git fetch` + fast-forward to `origin/domain_split`, leaving the submodule on that branch). No file edits, no local commits, no patches, not even temporary ones. When something in that pipeline needs to change, say what the change is and stop — the user makes it in `~/Programming/ppi-splitting-pipeline`, pushes, and then it is pulled here. A local edit is invisible upstream, is destroyed by the next pull, and its `bin/` is baked into the Docker image, so it produces a silently _wrong_ run rather than a failure. Bumping the pin shows up in this repo as `M subworkflows/external/ppi-splitting`, which looks like a content edit — say explicitly that only the pointer moved, and report the commit id.
- Nextflow reads only the root project's `nextflow.config`, so that pipeline's params come in via `includeConfig` of its `conf/params.config` and everything else it needs (container, per-process resources, the `gurobi` label, publishing) is restated in `conf/ppi_splitting.config`. Its `bin/` is baked into its Docker image, so the image tag and the submodule pin move together.
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
