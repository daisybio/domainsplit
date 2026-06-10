process BUILD_PPI_NEGATIVE_POOL {
    tag "build_ppi_negative_pool"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path negative_ppi_parquet
    val  min_n_tested
    val  self_interaction

    output:
    path "domainsplit.sqlite3",        emit: domainsplit_db
    path "neg_pool.npz",               emit: neg_pool
    path "uniprot_pfam_mapping.json",  emit: pfam_mapping
    path "versions.yml",               emit: versions

    script:
    def no_self = self_interaction ? "" : "--no-self"
    """
    cp "${domainsplit_db_in}" domainsplit.sqlite3

    build_ppi_negative_pool.py \\
        --db domainsplit.sqlite3 \\
        --parquet "${negative_ppi_parquet}" \\
        --pfam-mapping-out uniprot_pfam_mapping.json \\
        --pool-out neg_pool.npz \\
        --min-n-tested ${min_n_tested} \\
        ${no_self}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import sys; print(sys.version.split()[0])')
        pyarrow: \$(python3 -c 'import pyarrow; print(pyarrow.__version__)')
        numpy: \$(python3 -c 'import numpy; print(numpy.__version__)')
        sqlite3: \$(python3 -c 'import sqlite3; print(sqlite3.sqlite_version)')
    END_VERSIONS
    """
}
