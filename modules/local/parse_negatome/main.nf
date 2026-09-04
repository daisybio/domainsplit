process PARSE_NEGATOME {
    tag "parse_negatome"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path negatome_txt

    output:
    path "negatome_ddis.tsv", emit: ddis
    path "versions.yml",      emit: versions

    script:
    """
    parse_negatome.py \\
        --negatome ${negatome_txt} \\
        --out negatome_ddis.tsv \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
