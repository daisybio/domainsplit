process ANALYZE_DDI_BIAS {
    tag "bias_analysis"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path "domainsplit.sqlite3"

    output:
    path "bias_analysis", emit: report_dir
    path "versions.yml",  emit: versions

    script:
    """
    analyze_ddi_bias.py \
        --db domainsplit.sqlite3 \
        --outdir bias_analysis

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import sys; print(sys.version.split()[0])')
        numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
        matplotlib: \$(python3 -c 'import matplotlib; print(matplotlib.__version__)')
    END_VERSIONS
    """
}
