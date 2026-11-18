// reports/run_manifest.json -- what this run was made from and what it made.
// DATASET_AUDIT.md F10: the audited run could not be tied to its inputs, seed or
// per-split class ratios once its work directories were gone. See the script's
// docstring for what each field means.
process WRITE_RUN_MANIFEST {
    tag "run_manifest"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path base_json
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path hippie_tsv,           stageAs: 'manual/hippie/*'
    path ppidm_tsv,            stageAs: 'manual/ppidm/*'
    path negative_ppi_parquet, stageAs: 'manual/negative_ppi/*'
    // COMPUTE_PPIDM's ppidm_report.json and BUILD_PPIDM_SOURCES' ppidm_sources.tsv
    // (both also published under reports/); [] when PPIDM is not computed.
    path ppidm_reports,        stageAs: 'ppidm/*'
    // BUILD_CANDIDATE_NETWORK's gene_pfam_mapping.json; [] without `ilp_candidates`.
    path candidate_pfam_mapping, stageAs: 'candidate/*'

    output:
    path "run_manifest.json", emit: manifest
    path "versions.yml",      emit: versions

    script:
    def hippie_arg = hippie_tsv ? "--input hippie=${hippie_tsv}" : ''
    def ppidm_arg  = ppidm_tsv ? "--input ppidm_tsv=${ppidm_tsv}" : ''
    def neg_arg    = negative_ppi_parquet ? "--input negative_ppi_parquet=${negative_ppi_parquet}" : ''
    def ppidm_args = [ppidm_reports].flatten().findAll { f -> f }.collect { f -> "--input ${f.name}=${f}" }.join(' ')
    def mapping_arg = candidate_pfam_mapping ? "--input candidate_gene_pfam_mapping=${candidate_pfam_mapping}" : ''
    """
    write_run_manifest.py \\
        --base ${base_json} \\
        --db "${domainsplit_db_in}" \\
        ${hippie_arg} \\
        ${ppidm_arg} \\
        ${neg_arg} \\
        ${ppidm_args} \\
        ${mapping_arg} \\
        --out run_manifest.json \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
