// modules/local/structural/merge_scores/main.nf
process MERGE_SCORES {
    tag "merge_scores"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path dbstruct, stageAs: 'input.dbstruct.sqlite3'
    path scores_tsvs
    path confirmed_tsvs

    output:
    path "dbstruct.sqlite3", emit: dbstruct
    path "versions.yml",     emit: versions

    script:
    """
    merge_scores.py \\
        --db_in input.dbstruct.sqlite3 --db_out dbstruct.sqlite3 \\
        --scores ${scores_tsvs} --confirmed ${confirmed_tsvs} \\
        --versions versions.yml --process_name "${task.process}"
    """
}