// Included twice: as PARSE_PPIDM (positives) and as PARSE_INVERSE_PPIDM, whose
// `ext.args` (`--source-prefix INVERSE_PPIDM --negative`) and `ext.prefix`
// (`inverse_ppidm`) are set in conf/modules.config. The prefix keeps the two
// published TSVs from overwriting each other under external_ddis/.
process PARSE_PPIDM {
    tag "${task.ext.prefix ?: 'ppidm'}"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path ppidm_tsv
    val  classes

    output:
    path "*_ddis.tsv",     emit: ddis
    path "versions.yml",   emit: versions

    script:
    def prefix = task.ext.prefix ?: 'ppidm'
    def args   = task.ext.args ?: ''
    """
    parse_ppidm.py \\
        --ppidm ${ppidm_tsv} \\
        --classes "${classes}" \\
        ${args} \\
        --out ${prefix}_ddis.tsv \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
