process INSERT_PPI_NEGATIVE_SELECTION {
    tag "insert_ppi_negative_selection"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path neg_pool
    path pairs_deletion
    path pairs_random_addition
    path score_deletion
    path score_random_addition

    output:
    path "domainsplit.sqlite3",            emit: domainsplit_db
    path "negative_ppi_method_scores.tsv", emit: scores
    path "versions.yml",                   emit: versions

    script:
    """
    cp "${domainsplit_db_in}" domainsplit.sqlite3

    insert_ppi_negative_selection.py \\
        --db domainsplit.sqlite3 \\
        --pool "${neg_pool}" \\
        --pairs-deletion "${pairs_deletion}" \\
        --pairs-random-addition "${pairs_random_addition}" \\
        --score-deletion "${score_deletion}" \\
        --score-random-addition "${score_random_addition}" \\
        --scores-out negative_ppi_method_scores.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import sys; print(sys.version.split()[0])')
        numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
        sqlite3: \$(python3 -c 'import sqlite3; print(sqlite3.sqlite_version)')
    END_VERSIONS
    """

    stub:
    """
    touch domainsplit.sqlite3 negative_ppi_method_scores.tsv
    echo '"${task.process}":' > versions.yml
    echo '    stub: "true"' >> versions.yml
    """
}
