process ENRICH_STRUCTURAL_AF {
    tag "enrich_structural_af"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path db_struct, stageAs: 'input.dbstruct.sqlite3'
    path af_metadata
    val methods

    output:
    path "dbstruct.sqlite3", emit: dbstruct
    path "structures.h5", emit: structures
    path "structure_coverage.tsv", emit: coverage
    path "af3_to_model", emit: af3_to_model
    path "versions.yml", emit: versions

    script:
    """
    enrich_structural_af.py \\
        --db_in input.dbstruct.sqlite3 \\
        --db_out dbstruct.sqlite3 \\
        --af_metadata ${af_metadata} \\
        --structures_h5 structures.h5 \\
        --coverage structure_coverage.tsv \\
        --af3_to_model_dir af3_to_model \\
        --methods "${methods}" \\
        --max_mismatch_rate ${params.afdb_max_mismatch_rate} \\
        --versions versions.yml \\
        --process_name "${task.process}" \\
        --cpus ${task.cpus}
    """

    stub:
    """
    cp input.dbstruct.sqlite3 dbstruct.sqlite3
    touch structures.h5 structure_coverage.tsv
    mkdir af3_to_model
    touch af3_to_model/all_methods.tsv
    echo '"${task.process}"' > versions.yml
    echo '  stub: "true"' >> versions.yml
    """
}
