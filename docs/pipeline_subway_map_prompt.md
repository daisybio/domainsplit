# Claude Design prompt — `daisybio/domainsplit` subway map

## Prompt for Claude Design

Create a **subway-map-style overview graphic** for the `daisybio/domainsplit` Nextflow pipeline, in the visual idiom that nf-core uses for its pipeline metro maps (clean vector lines, rounded corners, station dots, monospaced station labels, light background, distinct line colors with a small legend).

### Conceptual model

- **Lines = split configurations** (not data sources). Every line traverses the same shared trunk of processing stations from left to right and only diverges near the right-hand end of the map, where each line picks one splitting strategy before all lines re-merge at the final `SPLIT_DATABASE` terminus.
- **Four lines** (in the legend, each its own color):
  1. **`random_ddi`** — random split on DDIs.
  2. **`random_denoise`** — random denoise split.
  3. **`random_discovery`** — random denoise split, discovery variant.
  4. **`minimal_leakage`** — leakage-minimised split. This line briefly **splits into two parallel branches** (protein and domain) over the divergence segment, then both branches re-merge at `SPLIT_DATABASE`.

### Layout

1. **Shared trunk (all 4 lines stacked / running parallel through the same stations).** Stations in order, left to right:
   - `INIT_DOMAINSPLIT_DB`
   - `DOWNLOAD_3DID_SQLITE`
   - `INSERT_DDIS` **and** `INSERT_NEGATIVE_DDIS` — short parallel pair (positive/negative DDI insertion), rejoining immediately.
   - `EXTRACT_UNIQUE_DOMAINS` — reads distinct Pfam IDs from the in-build DB's `domain_domain_interaction` table.
   - `DOWNLOAD_PFAM_ALIGNMENT`
   - `CREATE_PROTEIN_DOMAIN_MAPPING`
   - `DOWNLOAD_PROTT5_EMBEDDINGS`
   - `GENERATE_ESM_PROTEIN_EMBEDDINGS` **and** `GENERATE_ESM_DOMAIN_EMBEDDINGS` — render these as two stations on a short parallel pair (like an interchange) since the ESM subworkflow emits both. Both rejoin the trunk immediately after.
   - `INSERT_DOMAIN_GO_TERMS` → `INSERT_PROTEINS_WITH_EMBEDDINGS` → `INSERT_PROTEIN_GO_TERMS` → `INSERT_PPI` → `INSERT_DOMAIN_PROTEIN_MAPPING` — sequential `ENRICH_DDI_DATABASE` chain. Render as five close-spaced stations along the trunk; mark the last with a slightly larger dot to flag completion of `domainsplit.sqlite3`.
   - `EXTRACT_PROTEIN_SEQUENCES` **and** `EXTRACT_DOMAIN_SEQUENCES` — another short parallel pair off the trunk, rejoining immediately.
   - `MMSEQS_EASYCLUSTER` — interchange dot. Only the `minimal_leakage` line "uses" this station; the other three lines pass through without stopping (render as a hollow / pass-through dot for them, filled solid for `minimal_leakage`).

2. **Divergence segment** — the four lines fan out into their own split-method station, each clearly labelled with the line color:
   - `random_ddi` line → **`RANDOM_DDI_SPLIT`**
   - `random_denoise` line → **`RANDOM_DENOISE_SPLIT`** (denoise=false)
   - `random_discovery` line → **`RANDOM_DENOISE_SPLIT`** (denoise=true, discovery variant — label e.g. `RANDOM_DENOISE_SPLIT (discovery)`)
   - `minimal_leakage` line → branches into **`MINIMAL_LEAKAGE_SPLIT_PROTEIN`** and **`MINIMAL_LEAKAGE_SPLIT_DOMAIN`** (two stations on parallel sub-branches in the same line color).

3. **Re-merge + terminus** — all four lines converge into a single terminus station:
   - **`SPLIT_DATABASE`** — terminus, larger station marker, label the outputs as `split_databases/<method>/<split>.sqlite3` with `<split> ∈ {train, optimization, test}`.

### Visual notes

- Use a clean horizontal layout with one or two gentle bends; mimic the London Tube / nf-core map look.
- Station labels in a monospaced font, all uppercase to match Nextflow process names. Sub-labels (e.g. ratios, variant flags) in a smaller weight underneath.
- Legend in a corner: 4 lines with their color swatches and the split-method name.
- Add a small subtitle: "Pipeline: daisybio/domainsplit — domainsplit DB assembly + train/optimization/test splitting".
- Add a tiny annotation near `INSERT_DOMAIN_PROTEIN_MAPPING` ("completes domainsplit.sqlite3") and near `SPLIT_DATABASE` ("emits per-method, per-split SQLite DBs").

### What NOT to include

- No external input URLs as stations (3did, UniProt, Negatome, STRING, pfam2go are inputs to the `ENRICH_DDI_DATABASE` chain but should not appear on the map).
- No container, resource, or config detail.
- No samplesheet — this pipeline ignores it.

---

## Station notes (one short paragraph each)

### `DOWNLOAD_3DID_SQLITE`

**In:** 3did SQL dump (`url_3did`).
**Out:** `3did.sqlite3` — domain-domain interaction database in SQLite form.
**What:** Downloads the gzipped MySQL dump from 3did and converts it to SQLite using the vendored `bin/mysql2sqlite` awk script. This is the structural-DDI seed for the whole pipeline.

### `EXTRACT_UNIQUE_DOMAINS`

**In:** the in-build `domainsplit.sqlite3` (post-DDI-insertion).
**Out:** flat list of Pfam accession IDs (one per line).
**What:** Queries the in-build DB's `domain_domain_interaction` table for every Pfam ID present on either side of an inserted DDI (JOIN through `domain` table). Defines the universe of domains the pipeline will care about downstream.

### `DOWNLOAD_PFAM_ALIGNMENT`

**In:** Pfam ID list (fanned out one ID per task via `splitText`).
**Out:** one alignment file per Pfam ID.
**What:** Downloads the full alignment for each Pfam family from the InterPro Pfam API (`url_pfam_template`). Many parallel tasks, one per Pfam ID. Outputs are later `.collect()`-ed.

### `CREATE_PROTEIN_DOMAIN_MAPPING`

**In:** UniProt id-mapping table (`url_uniprot_id_mapping`) + all Pfam alignment files.
**Out:** protein <-> domain mapping table (CSV, `uniprot_id` keyed).
**What:** Joins UniProt accessions to Pfam domain occurrences using the alignments. Produces the master mapping of which proteins contain which Pfam domains — the backbone of the protein/domain bridge.

### `DOWNLOAD_PROTT5_EMBEDDINGS`

**In:** Distinct UniProt IDs derived from the protein↔domain mapping.
**Out:** ProtT5 per-residue embeddings (HDF5) for those proteins.
**What:** Fetches precomputed ProtT5 embeddings from the EBI mirror, restricted to the protein set the pipeline actually needs (rather than the full UniProt dump).

### `GENERATE_ESM_PROTEIN_EMBEDDINGS` / `GENERATE_ESM_DOMAIN_EMBEDDINGS`

**In:** UniProt SwissProt sequences (`url_uniprot_sequences`) + protein↔domain mapping.
**Out:** Two HDF5 files — one with per-protein ESM embeddings, one with per-domain ESM embeddings.
**What:** Runs the ESM3 model on a GPU to embed (a) full-length proteins and (b) the domain sub-sequences carved out by the Pfam mapping. Single Nextflow process, dual outputs. Requires `HF_TOKEN` for the gated EvolutionaryScale model.

### `ENRICH_DDI_DATABASE` (interchange / hub — sequential chain)

**In:** post-DDI `domainsplit.sqlite3` + pfam2go, UniProt sequences, protein↔domain mapping, ProtT5 embeddings, UniProt GO terms, STRING links, UniProt id-mapping, ESM protein embeddings, ESM domain embeddings.
**Out:** Single unified `domainsplit.sqlite3` (the enriched DDI database).
**What:** Five sequential `INSERT_*` processes, each opening the SQLite emitted by the previous step, committing one phase, and handing the DB forward — `INSERT_DOMAIN_GO_TERMS` → `INSERT_PROTEINS_WITH_EMBEDDINGS` → `INSERT_PROTEIN_GO_TERMS` → `INSERT_PPI` → `INSERT_DOMAIN_PROTEIN_MAPPING`. End state contains DDIs, PPIs (positive + negative), GO annotations, sequences, and both embedding modalities. This is the canonical artifact the splitting half of the pipeline consumes.

### `EXTRACT_PROTEIN_SEQUENCES` / `EXTRACT_DOMAIN_SEQUENCES`

**In:** `domainsplit.sqlite3`.
**Out:** Two FASTA files — `proteins.fasta` and `domains.fasta`.
**What:** Pulls the protein and domain sequences out of the domainsplit DB into FASTA so they can be clustered. Domain FASTAs are sub-sequences carved using the Pfam coordinates from `CREATE_PROTEIN_DOMAIN_MAPPING`.

### `MMSEQS_EASYCLUSTER` (interchange, used by `minimal_leakage` only)

**In:** Protein FASTA and domain FASTA.
**Out:** Two cluster TSVs — one for proteins, one for domains (`representative \t member`).
**What:** Clusters sequences at 40% identity / 80% coverage (`--min-seq-id 0.4 -c 0.8` via `ext.args` in `conf/modules.config`). Cluster assignments feed the `minimal_leakage` split methods so members of the same cluster can be forced into the same split (preventing sequence-similarity leakage). The other three split lines pass through this station without consuming its output.

### `RANDOM_DDI_SPLIT` — `random_ddi` line

**In:** `domainsplit.sqlite3` + split ratios `[train 0.6, optimization 0.2, test 0.2]`.
**Out:** Three TSVs of DDI IDs (one per split) + corresponding per-split path files.
**What:** Naive baseline split: shuffle DDIs, partition into train/optimization/test according to the ratios. No leakage control. Useful as a maximum-leakage upper bound.

### `RANDOM_DENOISE_SPLIT` — `random_denoise` line

**In:** `domainsplit.sqlite3`, denoise fraction `0.1`, discovery flag `false`.
**Out:** Per-split DDI ID files.
**What:** Random split with a denoising step — removes/holds out a fraction of low-confidence DDIs to clean the training signal. Discovery flag off, so kept DDIs stay in the training distribution.

### `RANDOM_DENOISE_SPLIT (discovery)` — `random_discovery` line

**In:** `domainsplit.sqlite3`, denoise fraction `0.1`, discovery flag `true`.
**Out:** Per-split DDI ID files.
**What:** Same module as `random_denoise` but `denoise=true` flips it into discovery mode: the held-out fraction is routed into test so the model is evaluated on its ability to recover putative novel interactions rather than memorise known ones.

### `MINIMAL_LEAKAGE_SPLIT_PROTEIN` / `MINIMAL_LEAKAGE_SPLIT_DOMAIN` — `minimal_leakage` line (two branches)

**In:** `domainsplit.sqlite3`, split ratios, and the corresponding cluster TSV (protein or domain) from MMseqs2.
**Out:** Per-split DDI ID files for the protein-leakage-minimised partition and the domain-leakage-minimised partition.
**What:** Simulated-annealing-style assignment that respects MMseqs2 cluster boundaries: every cluster is assigned to exactly one split so train/test cannot share near-identical sequences. The protein variant minimises protein-sequence leakage; the domain variant minimises domain-sequence leakage. Most expensive split.

### `SPLIT_DATABASE` (terminus)

**In:** All per-method, per-split DDI ID files (random_ddi, random_denoise, random_discovery, minimal_leakage_protein, minimal_leakage_domain) joined with their split labels and the full `domainsplit.sqlite3`.
**Out:** `split_databases/<method>/<split>.sqlite3` — one trimmed SQLite per (method × split) combination.
**What:** Materialisation step. For every (method, split) pair it produces a self-contained SQLite DB containing only the rows belonging to that split, ready to be consumed by a downstream training job. This is the pipeline's final published artifact.
