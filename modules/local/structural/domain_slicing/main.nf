process DOMAIN_SLICE {
    tag "domain_slice"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path pdb_dir_in, stageAs: 'input/pdb_files/'

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    """
    slice_domains.py \\
        --db_in  ${domainsplit_db_in} \\
        --db_out domainsplit.sqlite3 \\
        --pdb_dir ${pdb_dir_in}\\
        --versions versions.yml \\
        --process_name "${task.process}"
    """
}