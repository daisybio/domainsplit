// PPIDM (Alborzi et al. 2021, PLOS Comput Biol 17:e1008844) computed in-pipeline:
// domain pairs enriched in PPI co-occurrence (PPIDM) and, by the lower tail of
// the same hypergeometric test, pairs depleted in it (INVERSE_PPIDM). See
// PLAN_ppidm_in_pipeline.md for the method and its deviations from the paper.
//
// Every input is a local file staged by Nextflow from a `url_*` param (URL or
// path); none of these scripts downloads anything.

// One MITAB release file (BioGRID or IntAct zip) -> compact TSV. The IntAct
// file is 11.6 GB uncompressed and is streamed out of the zip, never extracted.
process EXTRACT_PPI_MITAB {
    tag "${format}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    tuple val(format), path(mitab)

    output:
    tuple val(format), path("${format}_ppis.tsv.gz"), emit: ppis
    path "versions.yml",                               emit: versions

    script:
    """
    extract_ppi_mitab.py \\
        --input ${mitab} \\
        --format ${format} \\
        --out ${format}_ppis.tsv.gz \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}

// The eight PPI sources, mapped onto the all-species reviewed UniProt -> Pfam
// map. That map is independent of `instance_tier` on purpose: the labels are the
// same in every tier, only BUILD_EXTERNAL_TEST's instance check differs.
process BUILD_PPIDM_SOURCES {
    tag "build_ppidm_sources"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path sprot_dats,      stageAs: 'sprot/*'
    path biogrid
    path intact
    path hippie,          stageAs: 'hippie/*'
    path string_detailed, stageAs: 'string/*'
    path sifts_chain_uniprot
    path did_flat

    output:
    path "source_*.tsv.gz",       emit: sources
    path "protein_pfam.tsv.gz",   emit: protein_pfam
    path "ppidm_sources.tsv",     emit: report
    path "versions.yml",          emit: versions

    script:
    def dat_args    = [sprot_dats].flatten().collect { d -> "--sprot-dat ${d}" }.join(' ')
    def string_args = [string_detailed].flatten().collect { s -> "--string-detailed ${s}" }.join(' ')
    """
    build_ppidm_sources.py \\
        ${dat_args} \\
        --biogrid ${biogrid} \\
        --intact ${intact} \\
        --hippie ${hippie} \\
        ${string_args} \\
        --sifts-chain-uniprot ${sifts_chain_uniprot} \\
        --did-flat ${did_flat} \\
        --hippie-min-score ${params.ppidm_hippie_min_score} \\
        --string-min-score ${params.ppidm_string_min_score} \\
        --intact-min-miscore ${params.ppidm_intact_min_miscore} \\
        --biogrid-min-pubs ${params.ppidm_biogrid_min_pubs} \\
        --outdir . \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}

// Consensus score, weights, threshold and both hypergeometric tails. The gold
// standard is the 3did-only database (after INSERT_3DID).
process COMPUTE_PPIDM {
    tag "compute_ppidm"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:2.0.0"

    input:
    path sources
    path protein_pfam
    path gold_db, stageAs: 'gold.domainsplit.sqlite3'

    output:
    path "ppidm_computed.tsv",         emit: ppidm
    path "inverse_ppidm_computed.tsv", emit: inverse_ppidm
    path "ppidm_scores.tsv.gz",        emit: scores
    path "ppidm_report.json",          emit: report
    path "versions.yml",               emit: versions

    script:
    def source_args = [sources].flatten().collect { s -> "--source ${s}" }.join(' ')
    def weights_arg = params.ppidm_source_weights ? "--weights '${params.ppidm_source_weights}'" : ''
    def thresh_arg  = params.ppidm_score_threshold != null ? "--threshold ${params.ppidm_score_threshold}" : ''
    """
    compute_ppidm.py \\
        --protein-pfam ${protein_pfam} \\
        ${source_args} \\
        --gold-db ${gold_db} \\
        --alpha ${params.ppidm_alpha} \\
        --weight-iterations ${params.ppidm_weight_iterations} \\
        ${weights_arg} \\
        ${thresh_arg} \\
        --seed ${params.seed} \\
        --outdir . \\
        --versions versions.yml \\
        --process-name "${task.process}"
    """
}
