process PREDICT_COMPLEX_RF {
    tag "predict_complex_rf"
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path meta_ppis

    output:
    path "pdb_files_rf/",          emit: pdb_files
    path "versions.yml",        emit: versions

    script:
    """
    predict_complex_rf.py \\
        --db_in  ${domainsplit_db_in} \\
        --outdir pdb_files_rf/ \\
        --fastatmp fasta_tmp \\
        --versions versions.yml \\
        --process_name "${task.process}"
    """
}
