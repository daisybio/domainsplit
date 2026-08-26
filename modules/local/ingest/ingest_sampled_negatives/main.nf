process INGEST_SAMPLED_NEGATIVES {
    tag "ingest_sampled_negatives"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    // `split_keys` are "<method>:<split>" strings aligned index-for-index with
    // `split_csvs`. The split identity travels in the channel rather than in the
    // file name: ppi-splitting only suffixes names when a row asks for more than
    // one negative set, and parsing that back would couple ingest to its naming
    // rule. Staging into numbered directories keeps two methods' `train.csv`
    // from colliding.
    val  split_keys
    path split_csvs, stageAs: '?/*'

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    // `.name` on a staged input already carries the `?/*` subdirectory
    // (e.g. `1/train.csv`), so it must not be prefixed with the index again.
    def split_args = [split_keys, split_csvs].transpose().collect { key, csv ->
        "--split ${key}:${csv.name}"
    }.join(' \\\n        ')
    """
    # NOT cp: on this cluster's NFS, coreutils uses copy_file_range(), which the
    # server satisfies as a server-side copy -- 567 MB "copied" in 0.8 s and then
    # the syscall never returns (BUILD_EXTERNAL_TEST hung 3 h with every byte
    # already on the server; --reflink=never does not opt out, coreutils still
    # takes that path). dd is a plain read()/write() loop, so the bytes actually
    # cross the wire and the call terminates.
    dd if="${domainsplit_db_in}" of=domainsplit.sqlite3 bs=4M status=none

    ingest_sampled_negatives.py \\
        --db domainsplit.sqlite3 \\
        ${split_args} \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
