process EXPORT_SPLIT_DDIS {
    tag "export_split_ddis"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path instances_tsv
    val  source

    output:
    path "split_ddis.csv", emit: ddis
    path "versions.yml",   emit: versions

    script:
    def drop_homomer = params.drop_homomer_ddi.toString().toBoolean() ? '--drop-homomer' : ''
    """
    export_split_inputs.py \\
        --mode ddis ${drop_homomer} \\
        --db "${domainsplit_db_in}" \\
        --instances ${instances_tsv} \\
        --source ${source} \\
        --out split_ddis.csv \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
