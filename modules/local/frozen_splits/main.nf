// Frozen splits (bin/frozen_splits.py).
//
// EXPORT_FROZEN_SPLITS  every normal run: what ppi-splitting handed to
//     INGEST_SPLITS plus FETCH_DOMAIN_META's outputs and the positive population,
//     as `frozen_splits/` with a SHA-256 manifest and the final master's
//     membership and protein fingerprints. Pass that directory as
//     `--frozen_splits` to reuse these exact splits in a later run.
// VERIFY_FROZEN_SPLITS  in a `--frozen_splits` run: bundle integrity, and this
//     run's rebuilt positive population must equal the frozen one. INGEST_SPLITS
//     waits on it.
// CHECK_FROZEN_DB  in a `--frozen_splits` run: the final master's membership
//     and protein fingerprints must equal the frozen run's. Structure and the
//     split databases wait on it.
// Both gates are `cache false`: a bundle file edited in place need not change
// what Nextflow hashes for the staged directory, so a cached pass proves nothing.

process EXPORT_FROZEN_SPLITS {
    tag "export_frozen_splits"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path fetch_files,    stageAs: 'fetch/*'        // sequences.fasta, go_annotations.tsv, species.tsv, instances.tsv, dropped_families.tsv
    path split_ddis,     stageAs: 'split_ddis/*'
    val  family_specs                               // [[method, split], ...] in the order of family_files
    path family_files,   stageAs: 'family??/*'
    val  instance_specs
    path instance_files, stageAs: 'instance??/*'
    val  params_json
    path master_db,      stageAs: 'master.sqlite3'  // the final master: fingerprints for CHECK_FROZEN_DB

    output:
    path "frozen_splits", emit: bundle
    path "versions.yml",  emit: versions

    script:
    def fetch_args    = [fetch_files].flatten().collect { f -> "--fetch ${f}" }.join(' ')
    def family_args   = [family_specs, [family_files].flatten()].transpose()
        .collect { ms, f -> "--family ${ms[0]}:${ms[1]}:${f}" }.join(' ')
    def instance_args = [instance_specs, [instance_files].flatten()].transpose()
        .collect { ms, f -> "--instance ${ms[0]}:${ms[1]}:${f}" }.join(' ')
    """
    cat > params.json <<'END_PARAMS'
    ${params_json}
    END_PARAMS
    frozen_splits.py export \\
        ${fetch_args} \\
        --split_ddis ${split_ddis} \\
        ${family_args} \\
        ${instance_args} \\
        --params "\$(cat params.json)" \\
        --db master.sqlite3 \\
        --out frozen_splits

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    mkdir frozen_splits
    echo '{}' > frozen_splits/bundle_manifest.json
    echo '"${task.process}":' > versions.yml
    echo '  stub: "true"' >> versions.yml
    """
}

// Integrity only, before anything reads the bundle: a tampered instances.tsv
// would otherwise shape (and publish) this run's split_ddis.csv first.
process VERIFY_FROZEN_BUNDLE {
    tag "verify_frozen_bundle"
    label 'process_single'
    cache false

    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path bundle, stageAs: 'bundle'

    output:
    stdout emit: ok

    script:
    """
    frozen_splits.py verify --bundle bundle
    """

    stub:
    """
    echo ok
    """
}

process VERIFY_FROZEN_SPLITS {
    tag "verify_frozen_splits"
    label 'process_single'
    cache false

    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path bundle, stageAs: 'bundle'
    path split_ddis, stageAs: 'split_ddis/*'

    output:
    stdout emit: ok

    script:
    """
    frozen_splits.py verify --bundle bundle --split_ddis ${split_ddis}
    """

    stub:
    """
    echo ok
    """
}

process CHECK_FROZEN_DB {
    tag "check_frozen_db"
    label 'process_single'
    cache false

    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path bundle, stageAs: 'bundle'
    path master_db, stageAs: 'master.sqlite3'

    output:
    stdout emit: ok

    script:
    """
    frozen_splits.py check_db --bundle bundle --db master.sqlite3
    """

    stub:
    """
    echo ok
    """
}
