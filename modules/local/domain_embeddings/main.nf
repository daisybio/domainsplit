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
    only the database knows `domain.id` -- see bin/export_domain_embeddings.py.
*/

include { SHARD_FASTA } from '../util/main.nf'

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
    container "docker.io/konstantinpelz/domainsplit-gpu:1.1.0"
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
    tuple val(model), path(input_fasta)

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
        ${require_gpu}
    """
}

// Collect one model's chunks and re-key them to {domain_id}/{instance_id}.
process EXPORT_DOMAIN_EMBEDDINGS {
    tag { model }
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

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
        --run-id "${workflow.sessionId}" \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}

workflow generate_domain_embeddings {
    take:
    domain_sequences   // ppi-splitting's sequences.fasta, keyed by instance id
    domainsplit_db     // the pruned master; supplies domain.id

    main:
    // Validated here rather than in the schema: nf-schema can check the string's
    // shape but not that a name has a `embedding_batch_size_<model>` param behind
    // it, and an unknown name would otherwise reach the task as a null batch size.
    def models = params.embedding_models_domainbench.toString().tokenize(',')*.trim().findAll { m -> m }
    def known = ['esm3', 'esmc', 'prott5']
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

        chunks = GENERATE_DOMAIN_EMBEDDINGS_CHUNK(channel.fromList(models).combine(shards))

        // No `size:` on the groupTuple: shard_fasta.py emits min(num_shards,
        // records) shards, so a fixture smaller than `embedding_shards` would hang
        // a sized group. Every chunk has to be finished before the export anyway.
        embeddings = EXPORT_DOMAIN_EMBEDDINGS(
            chunks.chunk.groupTuple().combine(domainsplit_db)
        ).embeddings

        ch_versions = Channel.empty().mix(
            SHARD_FASTA.out.versions,
            GENERATE_DOMAIN_EMBEDDINGS_CHUNK.out.versions,
            EXPORT_DOMAIN_EMBEDDINGS.out.versions,
        )
    }

    emit:
    embeddings
    versions = ch_versions
}
