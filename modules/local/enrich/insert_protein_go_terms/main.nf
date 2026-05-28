process INSERT_PROTEIN_GO_TERMS {
    tag "insert_protein_go_terms"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path uniprot_go_terms

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    """
    insert_protein_go_terms.py \\
        --db-in ${domainsplit_db_in} \\
        --uniprot-go-terms ${uniprot_go_terms} \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
