// Per-method domain_metadata.csv (bin/export_domain_metadata.py): ProtParam
// sequence features per ddi_split_membership instance. Needs no structure, so
// it reads the master before ENRICH_STRUCTURAL and stays cached across
// structure re-runs.
process EXPORT_DOMAIN_METADATA {
    tag "export_domain_metadata"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path db, stageAs: 'input.domainsplit.sqlite3'

    output:
    path "metadata/*/domain_metadata.csv",       emit: csv
    path "metadata/domain_metadata_report.tsv",  emit: report
    path "versions.yml",                         emit: versions

    script:
    """
    export_domain_metadata.py \\
        --db input.domainsplit.sqlite3 \\
        --outdir metadata \\
        --report metadata/domain_metadata_report.tsv \\
        --versions versions.yml \\
        --process_name "${task.process}"
    """

    stub:
    """
    mkdir -p metadata/stub
    touch metadata/stub/domain_metadata.csv metadata/domain_metadata_report.tsv
    echo '"${task.process}":' > versions.yml
    echo '  stub: "true"' >> versions.yml
    """
}
