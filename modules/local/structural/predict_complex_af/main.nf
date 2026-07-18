process PREDICT_COMPLEX_AF {
    tag "predict_complex_af"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path dbstruct, stageAs: 'input.dbstruct.sqlite3'

    output:
    path "pdb_files_af/",          emit: pdb_files
    path "versions.yml",        emit: versions

    script:
    """
    predict_complex_af.py \\
        --db_in  ${dbstruct} \\
        --outdir pdb_files_af/ \\
        --versions versions.yml \\
        --process_name "${task.process}"
    """
}
