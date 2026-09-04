process INSERT_3DID {
    tag "insert_3did"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path sqlite_3did

    output:
    path "domainsplit.sqlite3",   emit: domainsplit_db
    // The pairs 3did offered, before this step's dedup and before the
    // human-instance filter downstream. It is the only place that number exists
    // -- the DB keeps unique pairs -- and REPORT_DDI_ATTRITION's `offered`
    // column for 3did comes from here.
    path "insert_3did_counts.tsv", emit: counts
    path "versions.yml",           emit: versions

    script:
    """
    # NOT cp: on this cluster's NFS, coreutils uses copy_file_range(), which the
    # server satisfies as a server-side copy -- 567 MB "copied" in 0.8 s and then
    # the syscall never returns (BUILD_EXTERNAL_TEST hung 3 h with every byte
    # already on the server; --reflink=never does not opt out, coreutils still
    # takes that path). dd is a plain read()/write() loop, so the bytes actually
    # cross the wire and the call terminates.
    dd if="${domainsplit_db_in}" of=domainsplit.sqlite3 bs=4M status=none

    insert_3did.py \\
        --db domainsplit.sqlite3 \\
        --sqlite-3did ${sqlite_3did} \\
        --counts-out insert_3did_counts.tsv \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
