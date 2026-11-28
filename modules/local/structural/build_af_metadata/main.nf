// af_metadata.csv from the AFDB download manifest, AF3 job directories and/or a
// hand-made af_metadata CSV: deduplicated, validated, downloaded-before-inferred,
// with a stat fingerprint of every model and confidence file (bin/build_af_metadata.py).
//
// Two processes, because new AF3 jobs land in the scanned directories without
// any Nextflow input changing:
//   SCAN_AF_METADATA   `cache false`: builds the CSV on every run, emits its
//                      SHA-256 and publishes this run's inventory and counts
//                      (they change with unfinished or shadowed jobs even when
//                      the selected models, and so the digest, do not).
//   BUILD_AF_METADATA  cached on that digest: an unchanged scan reuses it, so the
//                      CSV keeps its work-dir path and ENRICH_STRUCTURAL_AF and
//                      everything after it stay cached. It fails if its own result
//                      differs from the scan (jobs added between the two).
// The job directories are read in place, never staged: they must be visible at
// the same path inside the container, like the model paths af_metadata holds.

def afMetadataArgs(manifest, af_metadata, inferred_dirs) {
    [
        manifest ? "--manifest ${manifest}" : '',
        af_metadata ? "--af_metadata ${af_metadata}" : '',
        inferred_dirs.collect { d -> "--inferred_dir '${d}'" }.join(' '),
    ].join(' ')
}

process SCAN_AF_METADATA {
    tag "scan_af_metadata"
    label 'process_single'
    cache false

    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path manifest,    stageAs: 'manifest/*'     // [] when not given
    path af_metadata, stageAs: 'af_metadata/*'  // [] when not given
    val  inferred_dirs                           // list of absolute directory paths, may be empty

    output:
    env 'AF_DIGEST',               emit: digest
    path "af_inventory.tsv",       emit: inventory
    path "af_metadata_report.tsv", emit: report

    script:
    """
    build_af_metadata.py ${afMetadataArgs(manifest, af_metadata, inferred_dirs)} \\
        --out af_metadata.csv \\
        --inventory af_inventory.tsv \\
        --report af_metadata_report.tsv \\
        --digest af_metadata_digest.txt
    AF_DIGEST=\$(cat af_metadata_digest.txt)
    """

    stub:
    """
    AF_DIGEST=stub
    touch af_inventory.tsv af_metadata_report.tsv
    """
}

process BUILD_AF_METADATA {
    tag "build_af_metadata"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path manifest,    stageAs: 'manifest/*'
    path af_metadata, stageAs: 'af_metadata/*'
    val  inferred_dirs
    val  expected_digest

    output:
    path "af_metadata.csv",        emit: af_metadata
    path "versions.yml",           emit: versions

    script:
    """
    build_af_metadata.py ${afMetadataArgs(manifest, af_metadata, inferred_dirs)} \\
        --out af_metadata.csv \\
        --inventory af_inventory.tsv \\
        --report af_metadata_report.tsv \\
        --digest af_metadata_digest.txt \\
        --expect_digest ${expected_digest}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    echo 'uniprot_id_a,uniprot_id_b,path,source' > af_metadata.csv
    echo '"${task.process}":' > versions.yml
    echo '  stub: "true"' >> versions.yml
    """
}
