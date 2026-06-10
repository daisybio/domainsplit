process SELECT_PPI_NEGATIVE_DANS {
    tag "select_ppi_negative_dans:seed=${seed}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    tuple val(seed), path(neg_pool)

    output:
    tuple val(seed), path("score_${seed}.json"), path("pairs_${seed}.tsv"), emit: result
    path "versions.yml", emit: versions

    script:
    """
    select_ppi_negative_dans.py \\
        --pool "${neg_pool}" \\
        --seed ${seed} \\
        --score-out score_${seed}.json \\
        --pairs-out pairs_${seed}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import sys; print(sys.version.split()[0])')
        numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
    END_VERSIONS
    """
}
