process PARSE_PPIDM {
    tag "parse_ppidm"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path ppidm_tsv
    val  classes

    output:
    path "ppidm_ddis.tsv", emit: ddis
    path "versions.yml",   emit: versions

    script:
    """
    parse_ppidm.py \\
        --ppidm ${ppidm_tsv} \\
        --classes "${classes}" \\
        --out ppidm_ddis.tsv \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
