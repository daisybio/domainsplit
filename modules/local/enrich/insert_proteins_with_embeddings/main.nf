process INSERT_PROTEINS_WITH_EMBEDDINGS {
    tag "insert_proteins_with_embeddings"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path uniprot_database
    path protein_domain_map
    path prott5_embeddings
    path esm_protein_embeddings

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    """
    insert_proteins_with_embeddings.py \\
        --db-in ${domainsplit_db_in} \\
        --uniprot-db ${uniprot_database} \\
        --protein-domain-map ${protein_domain_map} \\
        --prott5-embeddings ${prott5_embeddings} \\
        --esm-protein-embeddings ${esm_protein_embeddings} \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
