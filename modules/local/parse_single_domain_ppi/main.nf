process PARSE_SINGLE_DOMAIN_PPI {
    tag "parse_single_domain_ppi"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path hippie_tsv
    path swissprot_map
    val  min_score

    output:
    path "single_domain_ppi_ddis.tsv", emit: ddis
    path "versions.yml",               emit: versions

    script:
    """
    parse_single_domain_ppi.py \\
        --hippie ${hippie_tsv} \\
        --swissprot-map ${swissprot_map} \\
        --min-score ${min_score} \\
        --out single_domain_ppi_ddis.tsv \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
