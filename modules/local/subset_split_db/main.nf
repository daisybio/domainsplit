process SUBSET_SPLIT_DB {
    tag "${meta.method}_${meta.split}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    tuple val(meta), path(domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3')

    output:
    tuple val(meta), path("${meta.method}_${meta.split}.sqlite3"), emit: split_db
    path "versions.yml",                                          emit: versions

    script:
    def out_db = "${meta.method}_${meta.split}.sqlite3"
    """
    cp "${domainsplit_db_in}" ${out_db}

    subset_split_db.py \\
        --db ${out_db} \\
        --method ${meta.method} \\
        --split ${meta.split} \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
