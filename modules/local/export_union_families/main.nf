process EXPORT_UNION_FAMILIES {
    tag "export_union_families"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path external_ddis

    output:
    path "families.txt", emit: families
    path "versions.yml", emit: versions

    script:
    """
    export_split_inputs.py \\
        --mode families \\
        --db "${domainsplit_db_in}" \\
        --external-ddis ${external_ddis} \\
        --out families.txt \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
