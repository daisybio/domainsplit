// InterPro counters per membership Pfam family (bin/fetch_interpro_metadata.py),
// cached under <cache_dir>/interpro-<release>/. Needs network unless
// --interpro_release pins a non-live release whose cache is complete.
process FETCH_INTERPRO_METADATA {
    tag "interpro_${release ?: 'live'}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path db, stageAs: 'input.domainsplit.sqlite3'
    val cache_dir
    val release

    output:
    path "metadata/interpro_counts.tsv",    emit: tsv
    path "metadata/interpro_metadata.json", emit: summary
    path "versions.yml",                    emit: versions

    script:
    """
    mkdir -p metadata
    fetch_interpro_metadata.py \\
        --db input.domainsplit.sqlite3 \\
        --cache_dir "${cache_dir}" \\
        --release "${release}" \\
        --out metadata/interpro_counts.tsv \\
        --summary metadata/interpro_metadata.json \\
        --versions versions.yml \\
        --process_name "${task.process}"
    """

    stub:
    """
    mkdir -p metadata
    printf 'pfam_id\\n' > metadata/interpro_counts.tsv
    echo '{}' > metadata/interpro_metadata.json
    echo '"${task.process}":' > versions.yml
    echo '  stub: "true"' >> versions.yml
    """
}
