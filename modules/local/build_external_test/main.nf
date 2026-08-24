process BUILD_EXTERNAL_TEST {
    tag "build_external_test"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    val  methods            // method directories this test set is written under
    val  split_name

    output:
    path "domainsplit.sqlite3",           emit: domainsplit_db
    path "external_test_dropped.tsv",     emit: dropped
    path "versions.yml",                  emit: versions

    script:
    def method_args = methods.collect { method -> "--method ${method}" }.join(' ')
    """
    cp "${domainsplit_db_in}" domainsplit.sqlite3

    build_external_test.py \\
        --db domainsplit.sqlite3 \\
        ${method_args} \\
        --split ${split_name} \\
        --target ${params.ddi_examples_target} \\
        --pool-factor ${params.ddi_examples_pool_factor} \\
        --seed ${params.seed} \\
        --dropped-out external_test_dropped.tsv \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
