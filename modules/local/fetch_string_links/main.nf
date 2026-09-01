// Per-organism STRING links for every species the run's protein universe holds.
//
// Sits between PRUNE_UNREPRESENTED_DDIS and ENRICH_DDI_DATABASE: the `protein`
// table is final by then (ingest created the rows, the prune removed DDIs but no
// proteins), and INSERT_PPI is the consumer.
//
// Skipped entirely when `--url_string` names one explicit links file, which is
// how `-profile test` stays offline.
process FETCH_STRING_LINKS {
    tag "fetch_string_links"
    label 'process_low'
    // One HTTP request per taxon, so a transient failure is the common one; the
    // script already tolerates a 404 per taxon and only fails when *nothing*
    // downloaded, which is exactly the case worth retrying.
    label 'error_retry'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db, stageAs: 'input.domainsplit.sqlite3'
    path string_map
    // Absolute cache root, or '' for none. Not a `path`: the task *writes* to it,
    // and a staged relative path would resolve inside the work directory and be
    // deleted with it -- the same reason FETCH_DOMAIN_META takes it as a val.
    val cache_dir

    output:
    path "string_links.txt.gz",    emit: links
    path "string_taxa_report.tsv", emit: report
    path "versions.yml",           emit: versions

    script:
    def cache_arg = cache_dir ? "--cache-dir ${cache_dir}" : ''
    """
    fetch_string_links.py \\
        --db ${domainsplit_db} \\
        --string-map ${string_map} \\
        --out string_links.txt.gz \\
        --report string_taxa_report.tsv \\
        --release ${params.string_release} \\
        --min-proteins ${params.string_min_proteins_per_taxon} \\
        --jobs ${task.cpus} \\
        ${cache_arg} \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
