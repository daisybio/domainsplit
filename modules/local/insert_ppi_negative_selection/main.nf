process INSERT_PPI_NEGATIVE_SELECTION {
    tag "insert_ppi_negative_selection"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path score_jsons
    path pairs_tsvs

    output:
    path "domainsplit.sqlite3",          emit: domainsplit_db
    path "negative_ppi_seed_scores.tsv", emit: scores
    path "versions.yml",                 emit: versions

    script:
    """
    cp "${domainsplit_db_in}" domainsplit.sqlite3

    insert_ppi_negative_selection.py \\
        --db domainsplit.sqlite3 \\
        --scores-out negative_ppi_seed_scores.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import sys; print(sys.version.split()[0])')
        sqlite3: \$(python3 -c 'import sqlite3; print(sqlite3.sqlite_version)')
    END_VERSIONS
    """
}
