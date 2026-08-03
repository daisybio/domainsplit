process DOMAIN_SLICE {
    tag "domain_slice"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path dbstruct, stageAs: 'input.dbstruct.sqlite3'
    path pdb_dir_af, stageAs: 'input/pdb_files_af/'
    path pdb_dir_rf, stageAs: 'input/pdb_files_rf/'

    output:
    path "dbstruct.sqlite3", emit: dbstruct
    path "versions.yml",        emit: versions

    script:
    """
    slice_domains.py \\
        --db_in ${dbstruct} \\
        --db_out dbstruct.sqlite3 \\
        --pdb_dir_af ${pdb_dir_af}\\
        --pdb_dir_rf ${pdb_dir_rf}\\
        --versions versions.yml \\
        --process_name "${task.process}"
    """
}