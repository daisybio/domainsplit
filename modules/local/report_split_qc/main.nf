// Split QC on the final master (bin/report_split_qc.py): overlap between every
// method's splits, the ILP family-disjointness guarantee (fails the task when
// violated), 1:1 vs 1:10 test overlap, the external test set's train-degree
// profile, and -- with esm3_structure -- which membership instances have no
// structure embedding and why.
process REPORT_SPLIT_QC {
    tag "report_split_qc"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    val  ilp_methods       // comma-separated
    val  external_methods  // comma-separated
    path esm3_h5           // esm3_structure_domain_embeddings.h5, or []
    path afdb_report       // afdb_structures_report.tsv, or []

    output:
    // Named, not `*.tsv`: that would re-emit the staged afdb report.
    path "{split_overlap,family_leaks,test_variant_overlap,external_test_profile}.tsv", emit: reports
    path "esm3_structure_*.tsv", optional: true, emit: esm3_reports
    path "versions.yml", emit: versions

    script:
    // Read-only (mode=ro in the script), so used where it is staged.
    def esm3_args = esm3_h5 ? "--esm3-h5 ${esm3_h5} --afdb-report ${afdb_report}" : ''
    """
    report_split_qc.py \\
        --db "${domainsplit_db_in}" \\
        --ilp-methods "${ilp_methods}" \\
        --external-methods "${external_methods}" \\
        ${esm3_args} \\
        --outdir . \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
