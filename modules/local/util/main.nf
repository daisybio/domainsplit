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
process SHARD_FASTA {
    tag { "${input_fasta.simpleName}:${num_shards}" }
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

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

// Merge a collection of HDF5 chunks (plain or gzipped) into one HDF5 file.
// Nested groups are preserved via h5py's in-place copy.
//
// Caller passes:
//   - output_name : basename for the joined file (no extension)
//   - "chunk*"    : a `.collect()`-ed channel of chunk paths
process JOIN_HDF_FILES {
    tag { output_name }
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    val output_name
    path "chunk*"

    output:
    path "${output_name}.h5", emit: joined
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    from pathlib import Path
    import gzip
    import h5py
    import sys

    out_path = "${output_name}.h5"
    with h5py.File(out_path, "w") as out_h5:
        for chunk_file in sorted(Path().glob("chunk*")):
            print("Combining chunk file:", chunk_file)
            # accept both gzipped and plain HDF5 chunks
            try:
                opener = gzip.open(chunk_file, "rb")
                opener.peek(1)
                in_h5 = h5py.File(opener, "r")
            except (OSError, gzip.BadGzipFile):
                in_h5 = h5py.File(chunk_file, "r")
            with in_h5:
                for key in in_h5.keys():
                    in_h5.copy(key, out_h5)

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    h5py: {h5py.__version__}\\n")
    """
}
