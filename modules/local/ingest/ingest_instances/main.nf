process INGEST_INSTANCES {
    tag "ingest_instances"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path instances_tsv
    path sequences_fasta

    output:
    path "domainsplit.sqlite3",             emit: domainsplit_db
    path "protein_domain_mapping.csv.gz",   emit: protein_domain_map
    path "versions.yml",                    emit: versions

    script:
    """
    cp "${domainsplit_db_in}" domainsplit.sqlite3

    ingest_instances.py \\
        --db domainsplit.sqlite3 \\
        --instances ${instances_tsv} \\
        --sequences ${sequences_fasta} \\
        --mapping-out protein_domain_mapping.csv.gz \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
