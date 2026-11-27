// Per-method ddi_metadata.csv (bin/export_ddi_metadata.py): af_* interface
// scores per structured instance pair, scored from the slices and PAE blocks in
// structures.h5. Reads the master after STRUCTURE_COVERAGE_BALANCE: a pair
// masked in any split of a method is left out of that method's file. Without
// structure (structures_h5 = []) every file is header-only.
process EXPORT_DDI_METADATA {
    tag "export_ddi_metadata"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path db, stageAs: 'input.domainsplit.sqlite3'
    path structures_h5  // structures.h5, or []

    output:
    path "metadata/*/ddi_metadata.csv",         emit: csv
    path "metadata/ddi_metadata_coverage.tsv",  emit: coverage
    path "versions.yml",                        emit: versions

    script:
    def h5_arg = structures_h5 ? "--structures_h5 ${structures_h5}" : ''
    """
    export_ddi_metadata.py \\
        --db input.domainsplit.sqlite3 \\
        ${h5_arg} \\
        --pae_cutoff ${params.ddi_metadata_pae_cutoff} \\
        --contact_cutoff ${params.ddi_metadata_contact_cutoff} \\
        --outdir metadata \\
        --coverage metadata/ddi_metadata_coverage.tsv \\
        --versions versions.yml \\
        --process_name "${task.process}"
    """

    stub:
    """
    mkdir -p metadata/stub
    touch metadata/stub/ddi_metadata.csv metadata/ddi_metadata_coverage.tsv
    echo '"${task.process}":' > versions.yml
    echo '  stub: "true"' >> versions.yml
    """
}
