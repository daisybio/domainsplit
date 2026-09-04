#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Negative test for ppi-splitting's --split_only contract
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Under --split_only, PPI_SPLITTING requires every dataset's meta to say
    split_method = 'ilp' AND negative_sampling_method = 'ilp', and requires all
    five precomputed files. It validates and error()s rather than rewriting meta,
    because meta is the join key for every join()/combine(by: 0) in there.

    All of that enforcement used to live in buildDatasetsChannel(), which an
    embedding pipeline never calls -- so a bad embedded run produced a degenerate
    split in silence. That regression can only be caught from *this* side, which
    is why the test lives here and not in the submodule.

    Run (each case must abort with the dataset id in the message):

      nextflow run tests/contracts/split_only_bad_meta.nf \
          -c tests/contracts/nextflow.config --split_only true --ddi_mode true

      # control: the same meta without --split_only is not validated at all
      nextflow run tests/contracts/split_only_bad_meta.nf \
          -c tests/contracts/nextflow.config --ddi_mode true
      # -> no validation error; proceeds into CLUSTERING

      # control: drop a precomputed file -> the missing-files branch fires first
      nextflow run tests/contracts/split_only_bad_meta.nf \
          -c tests/contracts/nextflow.config --split_only true --ddi_mode true \
          --drop_partition true

    Nothing is executed: the error fires in a channel operator, before any task
    is submitted. `-preview` does NOT work -- it builds the DAG without running
    the dataflow, so the validation never sees an item. The staged files are
    placeholders for the same reason.
----------------------------------------------------------------------------*/

include { PPI_SPLITTING } from '../../subworkflows/external/ppi-splitting/main.nf'

// Deliberately mirrors metaOf() in workflows/domainsplit.nf, so that a change
// there which breaks the contract shows up here.
def badMeta() {
    return [
        id                       : 'external_test',
        embedding_model          : 'none',
        cdhit_identity           : params.cdhit_identity,
        cdhit_wordsize           : params.cdhit_wordsize,
        split_method             : params.bad_split_method,   // 'kahip', not 'ilp'
        edge_weight              : params.edge_weight,
        kahip_k                  : params.kahip_k,
        ilp_kahip_k              : params.ilp_kahip_k,
        train_split              : 0.9,
        val_split                : 0.1,
        test_split               : 0.0,
        ilp_epsilon              : params.ilp_epsilon,
        ilp_max_sec              : params.ilp_max_sec,
        negative_sampling_method : 'ilp',
        neg_ilp_time_limit       : params.neg_ilp_time_limit,
        neg_ilp_lambda_degree    : params.neg_ilp_lambda_degree,
        neg_ilp_lambda_taxon_pair: params.neg_ilp_lambda_taxon_pair,
        neg_ilp_lambda_self_loop : params.neg_ilp_lambda_self_loop,
        neg_ilp_lambda_jaccard   : 0,
    ]
}

// A function, not a closure local: closure locals of the workflow body are not
// visible from inside a nested closure.
def placeholder(name) {
    def f = file("${workDir}/split_only_contract/${name}")
    f.parent.mkdirs()
    f.text = "placeholder\n"
    return f
}

workflow {
    // An absent optional file is [], never null.
    def files = [
        ppis             : placeholder('ppis.csv'),
        sequences        : placeholder('sequences.fasta'),
        go_annotations   : [],
        species          : placeholder('species.tsv'),
        domain_instances : placeholder('instances.tsv'),
        blast_results    : [],
        candidate_network: [],
        partition        : params.drop_partition ? [] : placeholder('partition.txt'),
        node_mapping     : placeholder('node_mapping.txt'),
    ]

    PPI_SPLITTING( channel.of(tuple(badMeta(), files)) )
}
