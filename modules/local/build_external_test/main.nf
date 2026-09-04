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
    # The copy is the one step build_external_test.py cannot report on, and it is
    # where this task's first two runs stalled -- so bracket it. Those runs cloned
    # the *enriched* master; the workflow now schedules this before enrichment, so
    # the input is a few MB rather than ~590 MB of per-residue blobs.
    echo "[external_test] copy start \$(date -u +%FT%TZ) input=\$(stat -Lc %s "${domainsplit_db_in}") bytes"
    # NOT cp: on this cluster's NFS, coreutils uses copy_file_range(), which the
    # server satisfies as a server-side copy -- 567 MB "copied" in 0.8 s and then
    # the syscall never returns (BUILD_EXTERNAL_TEST hung 3 h with every byte
    # already on the server; --reflink=never does not opt out, coreutils still
    # takes that path). dd is a plain read()/write() loop, so the bytes actually
    # cross the wire and the call terminates.
    dd if="${domainsplit_db_in}" of=domainsplit.sqlite3 bs=4M status=none
    echo "[external_test] copy done  \$(date -u +%FT%TZ)"

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
