// The funnel ppi-splitting's own DDI attrition chart cannot show: its waterfall
// starts at the CSV EXPORT_SPLIT_DDIS handed it, so every DDI dropped *before*
// the splitting stage -- which under `instance_tier = human_only` is most of
// 3did -- is outside its frame. The counters existed but only in three
// `.command.log` files, which go away with the work directory.
//
// Reads files, not the DAG: every input here is already published, so the report
// can be regenerated from a finished run's `results/` without re-running
// anything. It takes the *pruned* master, so `surviving` means "in the published
// database".
process REPORT_DDI_ATTRITION {
    tag "report_ddi_attrition"
    label 'process_single'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path pruned_ddis
    path split_ddis
    path external_ddis
    path source_conflicts
    path offered_counts
    val  splitting_source

    output:
    path "ddi_source_attrition.tsv", emit: report
    path "versions.yml",             emit: versions

    script:
    // The DB is only ever read (mode=ro in the script), so it is used where it
    // is staged rather than copied -- unlike the modules that modify it, this
    // one has no reason to pay for a 500 MB dd.
    """
    report_ddi_attrition.py \\
        --db "${domainsplit_db_in}" \\
        --pruned-ddis ${pruned_ddis} \\
        --split-ddis ${split_ddis} \\
        --external-ddis ${external_ddis} \\
        --conflicts ${source_conflicts} \\
        --offered-counts ${offered_counts} \\
        --splitting-source ${splitting_source} \\
        --out ddi_source_attrition.tsv \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
