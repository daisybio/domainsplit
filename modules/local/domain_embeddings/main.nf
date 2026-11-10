/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Domain embeddings: sharded per-model inference over the domain FASTA, then
    one re-keyed HDF5 per model.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  union_sequences ── SHARD_FASTA ─┬─ GENERATE_DOMAIN_EMBEDDINGS_CHUNK (model × shard) ─┐
                                  └─ …                                                 │
                                      EXPORT_DOMAIN_EMBEDDINGS (per model) ────────────┘
                                          └─ <model>_domain_embeddings.h5

  - The input is ppi-splitting's own `sequences.fasta`, which is already the
    *domain* FASTA keyed by instance id (`fetch_domains.py` writes it under
    `{family}_{accession}_{start}_{end}`), so it is sharded directly. The old
    FILTER_SEQUENCES step re-derived both the domain FASTA and its keys from
    `protein_domain_mapping.csv.gz`; deleting it removes the whole class of
    key-mismatch bug rather than one instance of it.
  - Only the cut domain sequence is embedded, one pooled vector per instance, for
    every model in `params.embedding_models_domainbench`. Per-residue protein
    embeddings are retired: protein context correlates with the interaction
    partner, so a sliced protein embedding smuggles protein identity into a split
    that was partitioned on families.
  - One model per task, not all of them per task: independent retries and a
    per-model batch size, and a shard that OOMs costs one model.
  - EXPORT_DOMAIN_EMBEDDINGS takes the database as well as the chunks, because
    only the database maps an instance id to its Pfam family -- see
    bin/export_domain_embeddings.py. The published key is the Pfam accession,
    deliberately not `domain.id`: that surrogate differs between runs, so a file
    keyed on it read the wrong domain's vectors rather than failing.
*/

include { SHARD_FASTA } from '../util/main.nf'

// AlphaFold DB monomer coordinates for every parent protein of a domain
// instance, UniProt-numbered and checked against this run's UniProt sequences
// (bin/fetch_afdb_structures.py). CPU and network only: the GPU task that reads
// the H5 makes no network call. Runs only when `esm3_structure` is requested.
process FETCH_AFDB_STRUCTURES {
    tag "afdb_v${params.afdb_version}"
    label 'process_low'
    label 'error_retry'
    // The gpu env for ESM's own PDB parser (same atom37 layout ESM3 encodes);
    // nothing here touches a GPU, and no `process_gpu` label routes it to one.
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-gpu:2.0.0"

    input:
    path domain_sequences
    path uniprot_sequences
    // Absolute cache root or '' -- a val for the reason FETCH_STRING_LINKS gives.
    val cache_dir

    output:
    path "afdb_structures.h5",         emit: structures
    path "afdb_structures_report.tsv", emit: report
    path "versions.yml",               emit: versions

    script:
    def cache_arg = cache_dir ? "--cache-dir ${cache_dir}" : ''
    """
    fetch_afdb_structures.py \\
        --domain-fasta ${domain_sequences} \\
        --uniprot-fasta ${uniprot_sequences} \\
        --output-h5 afdb_structures.h5 \\
        --report afdb_structures_report.tsv \\
        --version ${params.afdb_version} \\
        --min-plddt ${params.afdb_min_plddt} \\
        --max-mismatch-rate ${params.afdb_max_mismatch_rate} \\
        --jobs ${task.cpus} \\
        ${cache_arg} \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}

// One pooled-embedding task per (model, shard).
process GENERATE_DOMAIN_EMBEDDINGS_CHUNK {
    tag { "${model}:${input_fasta.simpleName}" }
    // Both labels on purpose. `process_gpu` is what institutional configs key
    // their GPU queue and --gpus request on (nf-core/configs daisybio, for one);
    // `process_gpu_large` only carries our own cpu/memory/time sizing. Dropping
    // the first sent these tasks to the default CPU queue, where torch found no
    // GPU and a CPU ESM3 load got the task OOM-killed.
    label 'process_gpu'
    label 'process_gpu_large'
    secret 'HF_TOKEN'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-gpu:2.0.0"
    // Apptainer/Singularity `--env` requires KEY=VALUE (a bare name is rejected),
    // unlike Docker's `-e KEY` pass-through. HF_HOME / HUGGINGFACE_HUB_CACHE are
    // exported inside the task script from `params.embedding_hf_cache_dir`, so only
    // the HF_TOKEN secret has to cross the container boundary here. ProtT5 is
    // ungated and needs no token, but the process is one process for all models.
    containerOptions {
        workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer'
            ? '--env HF_TOKEN=$HF_TOKEN'
            : '-e HF_TOKEN'
    }

    input:
    // `structures` is the AFDB H5 for esm3_structure and `[]` for every other
    // model, so a new H5 never invalidates the cached esm3/esmc/prott5 tasks.
    tuple val(model), path(input_fasta), path(structures)

    output:
    tuple val(model), path("${input_fasta.simpleName}.${model}.h5"), emit: chunk
    path "versions.yml", emit: versions

    script:
    // Abort rather than fall back to CPU: a 1.4B model per task on CPU is not a
    // slow success, it is an OOM kill (exit 137) that the retry repeats.
    // `.toString().toBoolean()`: a CLI `--embedding_require_gpu false` arrives as the
    // String "false", which is truthy in Groovy, so the documented way to accept CPU
    // inference on a tiny input would otherwise not work.
    def require_gpu = params.embedding_require_gpu.toString().toBoolean() ? '--require-gpu' : ''
    def hf_cache    = params.embedding_hf_cache_dir ?: ''
    def batch_size  = params["embedding_batch_size_${model}"]
    def structure_args = model == 'esm3_structure'
        ? "--structures-h5 ${structures} --min-coverage ${params.esm3_structure_min_coverage}"
        : ''
    """
    if [ -n "${hf_cache}" ]; then
        mkdir -p "${hf_cache}"
        export HF_HOME="${hf_cache}"
        export HUGGINGFACE_HUB_CACHE="${hf_cache}"
    fi

    run_embeddings.py \\
        --input-fasta "${input_fasta}" \\
        --output-h5 "${input_fasta.simpleName}.${model}.h5" \\
        --versions versions.yml \\
        --process-name "${task.process}" \\
        --model ${model} \\
        --prott5-model "${params.embedding_prott5_model}" \\
        --batch-size ${batch_size} \\
        --max-len ${params.embedding_max_len} \\
        ${structure_args} \\
        ${require_gpu}
    """
}

// Collect one model's chunks and re-key them to {pfam_id}/{instance_id}.
process EXPORT_DOMAIN_EMBEDDINGS {
    tag { model }
    label 'process_medium'
    // CPU process sharing this module dir with the GPU one: general env, not ./environment.yml (gpu).
    conda "${moduleDir}/../../../docker/env/general.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    tuple val(model), path("chunk*"), path(domainsplit_db)

    output:
    path "${model}_domain_embeddings.h5", emit: embeddings
    path "versions.yml", emit: versions

    script:
    """
    export_domain_embeddings.py \\
        --db "${domainsplit_db}" \\
        --chunk-glob 'chunk*' \\
        --model ${model} \\
        --output-h5 "${model}_domain_embeddings.h5" \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}

workflow generate_domain_embeddings {
    take:
    domain_sequences   // ppi-splitting's sequences.fasta, keyed by instance id
    domainsplit_db     // the pruned master; supplies the instance -> Pfam family map
    uniprot_sequences  // PARSE_SWISSPROT's FASTA: the sequence AFDB models are checked against
    cache_dir          // absolute cache root or ''

    main:
    // Validated here rather than in the schema: nf-schema can check the string's
    // shape but not that a name has a `embedding_batch_size_<model>` param behind
    // it, and an unknown name would otherwise reach the task as a null batch size.
    def models = params.embedding_models_domainbench.toString().tokenize(',')*.trim().findAll { m -> m }
    def known = ['esm3', 'esmc', 'prott5', 'esm3_structure']
    def unknown = models.findAll { m -> !(m in known) }
    if (unknown) {
        error("--embedding_models_domainbench names unknown model(s) ${unknown.join(', ')} (known: ${known.join(', ')}).")
    }
    if (models.size() != new HashSet(models).size()) {
        error("--embedding_models_domainbench '${params.embedding_models_domainbench}' names a model more than once.")
    }

    // An empty model list is the skip switch, and it has to short-circuit before
    // SHARD_FASTA: with no models the fan-out is empty anyway, but sharding the
    // FASTA for nobody is a task and a staged copy for nothing.
    if (!models) {
        log.info "embedding_models_domainbench is empty -- no domain embeddings will be produced."
        embeddings = Channel.empty()
        ch_versions = Channel.empty()
    }
    else {
        shards = SHARD_FASTA(
            domain_sequences.map { fasta -> tuple([id: 'domain_sequences'], fasta) },
            params.embedding_shards,
        ).shards.flatten()

        // Fetched only when its one consumer is enabled.
        if ('esm3_structure' in models) {
            FETCH_AFDB_STRUCTURES(domain_sequences, uniprot_sequences, cache_dir)
            structures  = FETCH_AFDB_STRUCTURES.out.structures
            ch_versions = FETCH_AFDB_STRUCTURES.out.versions
        }
        else {
            // A placeholder, not `[]`: combine() would splice an empty list away.
            structures  = channel.value('none')
            ch_versions = Channel.empty()
        }

        chunks = GENERATE_DOMAIN_EMBEDDINGS_CHUNK(
            channel.fromList(models).combine(shards).combine(structures)
                .map { model, fasta, h5 -> tuple(model, fasta, model == 'esm3_structure' ? h5 : []) }
        )

        // No `size:` on the groupTuple: shard_fasta.py emits min(num_shards,
        // records) shards, so a fixture smaller than `embedding_shards` would hang
        // a sized group. Every chunk has to be finished before the export anyway.
        embeddings = EXPORT_DOMAIN_EMBEDDINGS(
            chunks.chunk.groupTuple().combine(domainsplit_db)
        ).embeddings

        ch_versions = ch_versions.mix(
            SHARD_FASTA.out.versions,
            GENERATE_DOMAIN_EMBEDDINGS_CHUNK.out.versions,
            EXPORT_DOMAIN_EMBEDDINGS.out.versions,
        )
    }

    emit:
    embeddings
    versions = ch_versions
}
