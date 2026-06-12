process BUILD_SWISSPROT_PFAM_MAP {
    tag "swissprot_map"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    val url

    output:
    path "swissprot_pfam_map.json", emit: map
    path "versions.yml",            emit: versions

    script:
    """
    build_swissprot_pfam_map.py \\
        --url "${url}" \\
        --out swissprot_pfam_map.json \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """

    stub:
    """
    touch swissprot_pfam_map.json
    echo '"${task.process}":' > versions.yml
    echo '    stub: "true"' >> versions.yml
    """
}
