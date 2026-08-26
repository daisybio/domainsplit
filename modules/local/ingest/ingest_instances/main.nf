process INGEST_INSTANCES {
    tag "ingest_instances"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path instances_tsv
    path sequences_fasta

    output:
    path "domainsplit.sqlite3",             emit: domainsplit_db
    path "protein_domain_mapping.csv.gz",   emit: protein_domain_map
    path "versions.yml",                    emit: versions

    script:
    """
    # NOT cp: on this cluster's NFS, coreutils uses copy_file_range(), which the
    # server satisfies as a server-side copy -- 567 MB "copied" in 0.8 s and then
    # the syscall never returns (BUILD_EXTERNAL_TEST hung 3 h with every byte
    # already on the server; --reflink=never does not opt out, coreutils still
    # takes that path). dd is a plain read()/write() loop, so the bytes actually
    # cross the wire and the call terminates.
    dd if="${domainsplit_db_in}" of=domainsplit.sqlite3 bs=4M status=none

    ingest_instances.py \\
        --db domainsplit.sqlite3 \\
        --instances ${instances_tsv} \\
        --sequences ${sequences_fasta} \\
        --mapping-out protein_domain_mapping.csv.gz \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
