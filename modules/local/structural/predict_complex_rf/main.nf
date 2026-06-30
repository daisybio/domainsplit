process PREDICT_COMPLEX_RF {
    tag "predict_complex_rf"
    label 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'

    output:
    path "pdb_files/",          emit: pdb_files
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    """
    predict_complex_rf.py \\
        --db_in  ${domainsplit_db_in} \\
        --db_out domainsplit.sqlite3 \\
        --outdir pdb_files/ \\
        --fastatmp fasta_tmp \\
        --versions versions.yml \\
        --process_name "${task.process}"
    """
}
