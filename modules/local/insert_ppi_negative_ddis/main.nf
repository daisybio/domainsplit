process INSERT_PPI_NEGATIVE_DDIS {
    tag "insert_ppi_negative_ddis"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path negative_ppi_parquet
    val  min_n_tested
    val  sampling_strategy
    val  self_interaction

    output:
    path "domainsplit.sqlite3",        emit: domainsplit_db
    path "uniprot_pfam_mapping.json",  emit: pfam_mapping
    path "versions.yml",               emit: versions

    script:
    def no_self = self_interaction ? "" : "--no-self"
    """
    cp "${domainsplit_db_in}" domainsplit.sqlite3

    build_ppi_negative_ddis.py \\
        --db domainsplit.sqlite3 \\
        --parquet "${negative_ppi_parquet}" \\
        --pfam-mapping-out uniprot_pfam_mapping.json \\
        --min-n-tested ${min_n_tested} \\
        --sampling-strategy "${sampling_strategy}" \\
        ${no_self}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import sys; print(sys.version.split()[0])')
        pyarrow: \$(python3 -c 'import pyarrow; print(pyarrow.__version__)')
        sqlite3: \$(python3 -c 'import sqlite3; print(sqlite3.sqlite_version)')
    END_VERSIONS
    """
}
