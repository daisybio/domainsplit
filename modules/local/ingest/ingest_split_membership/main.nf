process INGEST_SPLIT_MEMBERSHIP {
    tag "ingest_split_membership"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    // See INGEST_SAMPLED_NEGATIVES for why the (method, split) identity travels
    // beside the files rather than in their names.
    val  split_keys
    path split_csvs, stageAs: '?/*'

    output:
    path "domainsplit.sqlite3",        emit: domainsplit_db
    path "unresolved_split_rows.tsv",  emit: unresolved
    path "versions.yml",               emit: versions

    script:
    // `.name` on a staged input already carries the `?/*` subdirectory
    // (e.g. `1/train.csv`), so it must not be prefixed with the index again.
    def split_args = [split_keys, split_csvs].transpose().collect { key, csv ->
        "--split ${key}:${csv.name}"
    }.join(' \\\n        ')
    """
    cp --reflink=auto "${domainsplit_db_in}" domainsplit.sqlite3

    ingest_split_membership.py \\
        --db domainsplit.sqlite3 \\
        ${split_args} \\
        --unresolved-out unresolved_split_rows.tsv \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
