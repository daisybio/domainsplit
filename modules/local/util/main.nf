/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Shared utility processes reusable across modules.

    These are kept here so callers can `include as <alias>` to invoke them
    multiple times in one workflow (Nextflow DSL2 requires aliases for
    multi-invocations of an included process).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Split a (gzipped) FASTA into N contiguous gzipped shards.
// Each shard is emitted into the `shards` output channel; downstream fan-out
// with `.flatten()` to parallelize per-shard work across the cluster.
// Emits min(num_shards, records) shards, so a caller must not `groupTuple(size:)`
// on the fan-out with `num_shards` as the size.
process SHARD_FASTA {
    tag { "${input_fasta.simpleName}:${num_shards}" }
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    tuple val(meta), path(input_fasta)
    val num_shards

    output:
    path "${meta.id}_shard_*.fasta.gz", emit: shards
    path "versions.yml", emit: versions

    script:
    """
    shard_fasta.py \\
        --input-fasta "${input_fasta}" \\
        --output-prefix "${meta.id}_shard" \\
        --num-shards ${num_shards}

    python3 - <<'PY' > versions.yml
    import sys, Bio
    print('"${task.process}":')
    print(f"    python: {sys.version.split()[0]}")
    print(f"    biopython: {Bio.__version__}")
    PY
    """
}
