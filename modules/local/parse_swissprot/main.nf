process PARSE_SWISSPROT {
    tag "parse_swissprot"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path swissprot_dat
    val taxon_ids

    output:
    path "uniprot_go_terms.tsv.gz",    emit: go_terms
    path "swissprot_pfam.tsv.gz",      emit: pfam_tsv
    path "uniprot_sequences.fasta.gz", emit: sequences
    // Entry -> STRING id, from the entries' own `DR   STRING;` lines. Replaces the
    // per-organism <ORG>_<taxid>_idmapping.dat.gz download, and covers every
    // reviewed species the taxon filter admits rather than exactly one.
    path "uniprot_string_map.tsv.gz",  emit: string_map
    path "versions.yml",               emit: versions

    script:
    // An empty string keeps every species. Coerced through toString() because a
    // CLI `--swissprot_taxon_ids 9606` arrives as a String while the config
    // default may be a GString -- see the boolean-param note in CLAUDE.md.
    def taxa = taxon_ids.toString().trim()
    """
    parse_swissprot_dat.py \\
        --dat ${swissprot_dat} \\
        --go-terms uniprot_go_terms.tsv.gz \\
        --pfam-map swissprot_pfam.tsv.gz \\
        --sequences uniprot_sequences.fasta.gz \\
        --string-map uniprot_string_map.tsv.gz \\
        --taxon-ids "${taxa}" \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
