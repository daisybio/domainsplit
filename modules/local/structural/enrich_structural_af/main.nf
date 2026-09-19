process ENRICH_STRUCTURAL_AF {
    tag "enrich_structural_af"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path db_struct, stageAs: 'input.dbstruct.sqlite3'
    path af_metaddata

    output:
    path "dbstruct.sqlite3", emit: dbstruct
    path "structures.h5", emit: structures
    path "versions.yml", emit: versions

    script:
    """
    enrich_structural_af.py \\
        --db_in input.dbstruct.sqlite3 \\
        --db_out dbstruct.sqlite3 \\
        --af_metadata ${af_metaddata} \\
        --structures_h5 structures.h5 \\
        --versions versions.yml \\
        --process_name "${task.process}"
    """

    stub:
    """
    cp input.dbstruct.sqlite3 dbstruct.sqlite3
    touch structures.h5
    echo '"${task.process}"' > versions.yml
    echo '  stub: "true"' >> versions.yml
    """
}