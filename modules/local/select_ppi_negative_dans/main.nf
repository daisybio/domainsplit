process SELECT_PPI_NEGATIVE_DANS {
    tag "select_ppi_negative_dans:${method}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    val  method
    val  seed
    path neg_pool

    output:
    path "score_${method}.json", emit: score
    path "pairs_${method}.tsv",  emit: pairs
    path "versions.yml",         emit: versions

    script:
    """
    select_ppi_negative_dans.py \\
        --pool "${neg_pool}" \\
        --method ${method} \\
        --seed ${seed} \\
        --score-out score_${method}.json \\
        --pairs-out pairs_${method}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import sys; print(sys.version.split()[0])')
        numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
    END_VERSIONS
    """
}
