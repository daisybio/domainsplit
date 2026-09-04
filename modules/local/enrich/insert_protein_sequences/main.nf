process INSERT_PROTEIN_SEQUENCES {
    tag "insert_protein_sequences"
    // process_medium, not process_high: the memory this step used to need was the
    // per-residue HDF5 reads, and those are gone.
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path uniprot_database
    path protein_domain_map

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    """
    insert_protein_sequences.py \\
        --db-in ${domainsplit_db_in} \\
        --uniprot-db ${uniprot_database} \\
        --protein-domain-map ${protein_domain_map} \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
