#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    daisybio/domainsplit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/daisybio/domainsplit
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TYPED PARAMETER DECLARATIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Types only. The *defaults* stay where they already are -- this pipeline's own
    in nextflow.config, the imported ppi-splitting ones in their params.config --
    so nothing here duplicates a value and the two pipelines still cannot drift.

    Nextflow's v2 (strict) parser, the default since 26.x, stopped inferring types
    for command-line parameters: `--seed 7` arrives as the String "7" and
    `--ddi_mode false` as the String "false". Two different failures follow:

      * A param this pipeline's schema declares is rejected outright --
        "Value is [string] but should be [integer]".
      * A param in `validation.ignoreParams` (every imported ppi-splitting one)
        is not validated at all, so the String flows silently into the code.
        `--split_only false` then reads as truthy, which is why
        workflows/domainsplit.nf has to write
        `params.ppi_splitting_multi_negset.toString().toBoolean()`.

    Declaring the type makes Nextflow coerce the CLI value before either can
    happen. `validation.lenientMode` does NOT do this -- it only widens the other
    direction, letting an integer satisfy a string type -- and
    `NXF_SYNTAX_PARSER=v1` only works by reverting to the retired parser.

    Only the *entry* script's declarations count: a `params { }` block inside an
    included script (the ppi-splitting submodule's own main.nf) is ignored, which
    is why the imported params have to be repeated here. Verified, not assumed.

    `Float`, not `Double`/`BigDecimal`/`Number`: those reject the String outright
    instead of coercing it. Booleans are `Boolean`.

    The nf-core boilerplate flags (--help, --version, --monochrome_logs, ...) are
    deliberately absent: they are used as bare flags, and `help` is a
    string-or-boolean union that a type here would break.
----------------------------------------------------------------------------------------
*/
params {
    // ---- this pipeline's own (declared in nextflow_schema.json) ----
    seed: Integer
    negative_ppi_min_n_tested: Integer
    ddi_examples_target: Integer
    ddi_examples_pool_factor: Integer
    embedding_shards: Integer
    embedding_max_len: Integer
    embedding_batch_size_esm3: Integer
    embedding_batch_size_esmc: Integer
    embedding_batch_size_prott5: Integer
    hippie_min_score: Float
    split_train_fraction: Float
    split_val_fraction: Float
    split_test_fraction: Float
    ddi_mode: Boolean
    embedding_require_gpu: Boolean
    ppi_splitting_multi_negset: Boolean

    // ---- imported from ppi-splitting via
    //      `includeConfig 'subworkflows/external/ppi-splitting/conf/params.config'`.
    //      In validation.ignoreParams, so nothing else type-checks them.
    //      Keep in step with that list when the submodule is bumped.
    split_only: Boolean
    cdhit_identity: Float
    cdhit_wordsize: Integer
    heatmap_max_per_split: Integer
    kahip_seed: Integer
    kahip_k: Integer
    ilp_kahip_k: Integer
    train_split: Float
    val_split: Float
    test_split: Float
    ilp_epsilon: Float
    ilp_max_sec: Integer
    neg_ilp_lambda_degree: Float
    neg_ilp_lambda_taxon_pair: Float
    neg_ilp_lambda_self_loop: Float
    neg_ilp_lambda_jaccard: Float
    neg_ilp_time_limit: Integer
    neg_ilp_mip_gap: Float
    ddi_select_max_sec: Integer
    ddi_max_ilp_candidates: Integer
    ddi_lambda_diversity: Float
    ddi_shortlist_factor: Integer
    ddi_select_verbose: Boolean
    ddi_candidate_factor: Integer
}

include { DOMAINSPLIT  } from './workflows/domainsplit'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_domainsplit_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_domainsplit_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow DAISYBIO_DOMAINSPLIT {

    main:

    //
    // WORKFLOW: Run pipeline
    //
    DOMAINSPLIT ()

emit:
    domainsplit_db = DOMAINSPLIT.out.domainsplit_db
    split_db   = DOMAINSPLIT.out.split_db
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    DAISYBIO_DOMAINSPLIT ()

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
    )

    publish:
    domainsplit_db = DAISYBIO_DOMAINSPLIT.out.domainsplit_db
    split_db   = DAISYBIO_DOMAINSPLIT.out.split_db
}

output {
    domainsplit_db {
    }
    split_db {
        path {
            it[1] >> "databases/${it[0].method}/${it[0].split}.sqlite3"
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
