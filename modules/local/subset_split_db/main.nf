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
    # No copy step: subset_split_db.py creates ${out_db} and pulls only the
    # surviving rows out of the master, which it opens read-only. Cloning the
    # master first meant reading and writing ~99% per-residue embedding blobs
    # this split does not keep, once per (method, split).
    subset_split_db.py \\
        --source "${domainsplit_db_in}" \\
        --out ${out_db} \\
        --method ${meta.method} \\
        --split ${meta.split} \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
