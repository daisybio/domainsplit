process INSERT_DOMAIN_PROTEIN_MAPPING {
    tag "insert_domain_protein_mapping"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path protein_domain_map
    path esm_domain_embeddings

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    """
    insert_domain_protein_mapping.py \\
        --db-in ${domainsplit_db_in} \\
        --protein-domain-map ${protein_domain_map} \\
        --esm-domain-embeddings ${esm_domain_embeddings} \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
