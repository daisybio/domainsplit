// AF3 prediction queue (bin/build_af3_queue.py): every protein pair of the
// structure methods' instance pairs that has no AF model yet, in fixed level
// order, plus the pairs whose model could not be sliced at all. Reads the master
// straight after ENRICH_STRUCTURAL_AF, so the coverage-balance guard's masking
// never hides a model that exists.
process BUILD_AF3_QUEUE {
    tag "build_af3_queue"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path db_struct, stageAs: 'input.dbstruct.sqlite3'
    path af_metadata, stageAs: 'af_metadata/*'
    path af3_to_model, stageAs: 'gap/*'  // ENRICH_STRUCTURAL_AF's af3_to_model dir
    val methods

    output:
    // The whole published af3_to_model/: ENRICH_STRUCTURAL_AF's gap reports
    // copied in plus the queue files, so one task owns the directory.
    path "af3_to_model",   emit: dir
    path "af3_queue.json", emit: counts   // for run_manifest.json only
    path "versions.yml",   emit: versions

    script:
    """
    mkdir af3_to_model
    cp -L ${af3_to_model}/* af3_to_model/
    build_af3_queue.py \\
        --db input.dbstruct.sqlite3 \\
        --af_metadata ${af_metadata} \\
        --gap_report ${af3_to_model}/all_methods.tsv \\
        --methods "${methods}" \\
        --seed ${params.seed} \\
        --out_dir af3_to_model \\
        --counts af3_queue.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    mkdir af3_to_model
    cp -L ${af3_to_model}/* af3_to_model/
    touch af3_to_model/queue.tsv af3_to_model/queue_summary.txt af3_to_model/unusable_models.tsv
    echo '{}' > af3_queue.json
    echo '"${task.process}":' > versions.yml
    echo '  stub: "true"' >> versions.yml
    """
}
