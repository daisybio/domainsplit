// Coverage-balance guard (bin/balance_structure_coverage.py), after ANNOTATE_DDI:
// per (method, split), structure coverage of the two labels is equalised by
// masking has_structure / z_score, or the whole split's structures are masked
// when coverage is below structure_min_coverage. Off with
// --structure_balance_coverage false: nothing masked, report still written.
process STRUCTURE_COVERAGE_BALANCE {
    tag "structure_coverage_balance"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path db_struct, stageAs: 'input.dbstruct.sqlite3'

    output:
    path "dbstruct.sqlite3",                emit: dbstruct
    path "structure_coverage_balance.tsv",  emit: report
    path "structure_coverage_balance.json", emit: summary
    path "versions.yml",                    emit: versions

    script:
    """
    balance_structure_coverage.py \\
        --db_in input.dbstruct.sqlite3 \\
        --db_out dbstruct.sqlite3 \\
        --seed ${params.seed} \\
        --enabled ${params.structure_balance_coverage.toString().toBoolean()} \\
        --min_coverage ${params.structure_min_coverage} \\
        --report structure_coverage_balance.tsv \\
        --summary structure_coverage_balance.json \\
        --versions versions.yml \\
        --process_name "${task.process}"
    """

    stub:
    """
    cp input.dbstruct.sqlite3 dbstruct.sqlite3
    touch structure_coverage_balance.tsv
    echo '{}' > structure_coverage_balance.json
    echo '"${task.process}":' > versions.yml
    echo '  stub: "true"' >> versions.yml
    """
}
