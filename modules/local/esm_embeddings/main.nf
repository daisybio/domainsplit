/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ESM embeddings: FILTER_SEQUENCES + sharded ESM3 / ESMC inference + join.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  protein fasta ── SHARD_PROTEIN_FASTA ─┬─ GENERATE_PROTEIN_ESM_EMBEDDINGS_CHUNK ─┐
                                        └─ …                          (N tasks) └─ JOIN_PROTEIN_EMBEDDINGS ── esm_protein_embeddings.h5
  domain fasta  ── SHARD_DOMAIN_FASTA  ─┬─ GENERATE_DOMAIN_ESM_EMBEDDINGS_CHUNK ──┐
                                        └─ …                          (M tasks) └─ JOIN_DOMAIN_EMBEDDINGS  ── esm_domain_embeddings.h5

  - Chunk processes call bin/run_esm_embeddings.py (batched bf16 inference).
  - Protein chunks: per-residue (L+2, D) fp16 datasets (BOS + seq + EOS).
  - Domain chunks: GPU mean-pool across all tokens (BOS+seq+EOS) -> (D,) fp16.
*/

include { SHARD_FASTA as SHARD_PROTEIN_FASTA       } from '../util/main.nf'
include { SHARD_FASTA as SHARD_DOMAIN_FASTA        } from '../util/main.nf'
include { JOIN_HDF_FILES as JOIN_PROTEIN_EMBEDDINGS } from '../util/main.nf'
include { JOIN_HDF_FILES as JOIN_DOMAIN_EMBEDDINGS  } from '../util/main.nf'

process FILTER_SEQUENCES {
    tag { "${protein_domain_map.simpleName}" }
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path protein_domain_map
    path uniprotkb_database

    output:
    tuple val(protein_meta), path("uniprot_filtered.fasta.gz"), emit: protein_sequences
    tuple val(domain_meta), path("domain_sequences.fasta.gz"), emit: domain_sequences
    path "versions.yml", emit: versions

    script:
    protein_meta = [id: "protein_sequences"]
    domain_meta = [id: "domain_sequences"]
    """
    #!/usr/bin/env python3
    import csv
    import gzip
    from Bio import SeqIO

    uniprot_ids = set()
    with gzip.open("${protein_domain_map}", "rt") as protein_domain_map_text, \\
         gzip.open("domain_sequences.fasta.gz", "wt") as domain_sequences_fasta:
        reader = csv.DictReader(protein_domain_map_text)
        for row in reader:
            domain_sequences_fasta.write(
                f">{row['pfam_id']}_{row['uniprot_id']}_{row['start_pos']}_{row['end_pos']}\\n{row['sequence']}\\n"
            )
            uniprot_ids.add(row['uniprot_id'])

    with gzip.open("uniprot_filtered.fasta.gz", "wt") as uniprot_filtered_fasta, \\
         gzip.open("${uniprotkb_database}", "rt") as uniprotkb_fasta:
        for record in SeqIO.parse(uniprotkb_fasta, "fasta"):
            record.id = record.id.split("|")[1]
            if record.id in uniprot_ids:
                record.description = ""
                SeqIO.write(record, uniprot_filtered_fasta, "fasta")

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        import sys, Bio
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    biopython: {Bio.__version__}\\n")
    """
}

// Per-residue protein embeddings. One task per FASTA shard.
process GENERATE_PROTEIN_ESM_EMBEDDINGS_CHUNK {
    tag { input_fasta.simpleName }
    label 'process_gpu_large'
    secret 'HF_TOKEN'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-gpu:1.0.0"
    containerOptions {
        workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer'
            ? '--env HF_TOKEN --env HF_HOME --env HUGGINGFACE_HUB_CACHE'
            : '-e HF_TOKEN -e HF_HOME -e HUGGINGFACE_HUB_CACHE'
    }

    input:
    path input_fasta

    output:
    path "${input_fasta.simpleName}.esm.h5", emit: chunk
    path "versions.yml", emit: versions

    script:
    def smoke = params.esm_smoke_test ? 100 : 0
    def hf_cache = params.esm_hf_cache_dir ?: ''
    """
    if [ -n "${hf_cache}" ]; then
        mkdir -p "${hf_cache}"
        export HF_HOME="${hf_cache}"
        export HUGGINGFACE_HUB_CACHE="${hf_cache}"
    fi

    run_esm_embeddings.py \\
        --input-fasta "${input_fasta}" \\
        --output-h5 "${input_fasta.simpleName}.esm.h5" \\
        --versions versions.yml \\
        --process-name "${task.process}" \\
        --mode per_residue \\
        --batch-size ${params.esm_batch_size_protein} \\
        --max-len ${params.esm_max_len} \\
        --smoke-limit ${smoke}
    """
}

// GPU-pooled domain embeddings. One task per FASTA shard.
process GENERATE_DOMAIN_ESM_EMBEDDINGS_CHUNK {
    tag { input_fasta.simpleName }
    label 'process_gpu_large'
    secret 'HF_TOKEN'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-gpu:1.0.0"
    containerOptions {
        workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer'
            ? '--env HF_TOKEN --env HF_HOME --env HUGGINGFACE_HUB_CACHE'
            : '-e HF_TOKEN -e HF_HOME -e HUGGINGFACE_HUB_CACHE'
    }

    input:
    path input_fasta

    output:
    path "${input_fasta.simpleName}.esm.h5", emit: chunk
    path "versions.yml", emit: versions

    script:
    def smoke = params.esm_smoke_test ? 100 : 0
    def hf_cache = params.esm_hf_cache_dir ?: ''
    """
    if [ -n "${hf_cache}" ]; then
        mkdir -p "${hf_cache}"
        export HF_HOME="${hf_cache}"
        export HUGGINGFACE_HUB_CACHE="${hf_cache}"
    fi

    run_esm_embeddings.py \\
        --input-fasta "${input_fasta}" \\
        --output-h5 "${input_fasta.simpleName}.esm.h5" \\
        --versions versions.yml \\
        --process-name "${task.process}" \\
        --mode pooled \\
        --batch-size ${params.esm_batch_size_domain} \\
        --max-len ${params.esm_max_len} \\
        --smoke-limit ${smoke}
    """
}

workflow generate_esm_embeddings {
    take:
    uniprotkb_database
    protein_domain_map

    main:
    filter_result = FILTER_SEQUENCES(protein_domain_map, uniprotkb_database)

    protein_shards = SHARD_PROTEIN_FASTA(filter_result.protein_sequences, params.esm_protein_shards).shards.flatten()
    domain_shards  = SHARD_DOMAIN_FASTA(filter_result.domain_sequences,  params.esm_domain_shards ).shards.flatten()

    protein_chunks = GENERATE_PROTEIN_ESM_EMBEDDINGS_CHUNK(protein_shards)
    domain_chunks  = GENERATE_DOMAIN_ESM_EMBEDDINGS_CHUNK(domain_shards)

    protein_embeddings = JOIN_PROTEIN_EMBEDDINGS('esm_protein_embeddings', protein_chunks.chunk.collect()).joined
    domain_embeddings  = JOIN_DOMAIN_EMBEDDINGS('esm_domain_embeddings',  domain_chunks.chunk.collect() ).joined

    emit:
    protein_embeddings
    domain_embeddings
}
