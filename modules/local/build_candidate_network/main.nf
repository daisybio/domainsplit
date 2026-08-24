process BUILD_CANDIDATE_NETWORK {
    tag "build_candidate_network"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db
    path negative_ppi_parquet
    val  min_n_tested
    path gene_mapping    // a previous run's gene_pfam_mapping.json, or [] to query UniProt

    output:
    path "candidate_network.csv",    emit: candidate_network
    path "gene_pfam_mapping.json",   emit: pfam_mapping
    path "versions.yml",             emit: versions

    script:
    // The one live API call in this pipeline. Supplying the mapping removes it,
    // which is what lets `-profile test` run offline.
    def mapping_arg = gene_mapping ? "--mapping-in ${gene_mapping}" : ''
    """
    build_candidate_network.py \\
        --db "${domainsplit_db}" \\
        --parquet "${negative_ppi_parquet}" \\
        --mapping-out gene_pfam_mapping.json \\
        ${mapping_arg} \\
        --network-out candidate_network.csv \\
        --min-n-tested ${min_n_tested}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 -c 'import sys; print(sys.version.split()[0])')
        pyarrow: \$(python3 -c 'import pyarrow; print(pyarrow.__version__)')
        sqlite3: \$(python3 -c 'import sqlite3; print(sqlite3.sqlite_version)')
    END_VERSIONS
    """
}
