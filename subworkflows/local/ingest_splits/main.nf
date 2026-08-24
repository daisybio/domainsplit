/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INGEST_SPLITS -- pull ppi-splitting's results into the domainsplit SQLite.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Three steps, and the order between them is load-bearing:

      1. INGEST_INSTANCES        instances.tsv + sequences.fasta -> domain,
                                 protein, domain_protein_map (instance level).
                                 Also emits the protein_domain_mapping.csv.gz
                                 view that the ESM and ENRICH steps consume,
                                 replacing CURATE_DOMAINS' output byte-for-byte
                                 in shape.
      2. INGEST_SAMPLED_NEGATIVES  every label=0 family pair, from every split of
                                 every negative set, as source='sampled_negative'.
                                 Must run before INSERT_EXTERNAL_SOURCES -- see
                                 the insertion-order rule in bin/ddi_db_utils.py.
      3. INGEST_SPLIT_MEMBERSHIP records which DDIs and which instance pairs
                                 belong to which (method, split). Needs both
                                 previous steps: label=1 pairs resolve against
                                 3did, label=0 pairs against step 2, and both
                                 instances against step 1.

    Only four of ppi-splitting's outputs are read: instances.tsv,
    sequences.fasta, the family-level split CSVs and the instance-level
    `*_instances.csv`. `examples/*_universe.txt`, `_reserve.txt` and
    `unclaimed.txt` are ignored -- they are ppi-splitting's own bookkeeping.

    `val` -> `validation` renaming happens inside the ingest scripts
    (bin/split_io.py), at this boundary; ppi-splitting keeps `val` internally.

    Both split channels carry `[method, split, csv]`, where `method` is the
    output directory name -- the negative-set -> method mapping (`ilp` ->
    `minimal_leakage`, `ilp_candidates` -> `minimal_leakage_hcni`, ...) is the
    caller's job, because only the caller knows which dataset row produced which
    negative set.
----------------------------------------------------------------------------*/

include { INGEST_INSTANCES          } from '../../../modules/local/ingest/ingest_instances/main.nf'
include { INGEST_SAMPLED_NEGATIVES  } from '../../../modules/local/ingest/ingest_sampled_negatives/main.nf'
include { INGEST_SPLIT_MEMBERSHIP   } from '../../../modules/local/ingest/ingest_split_membership/main.nf'

// Freeze a `[method, split, file]` channel into a deterministically ordered pair
// of lists. `collect()` alone would gather in task-completion order, which would
// make the rendered command line -- and therefore the task hash -- unstable.
def orderedSplits(splits_ch) {
    return splits_ch.toSortedList { a, b -> "${a[0]}:${a[1]}" <=> "${b[0]}:${b[1]}" }
}

workflow INGEST_SPLITS {
    take:
    domainsplit_db_in
    instances_tsv       // ppi-splitting's instances.tsv, for the whole run
    sequences_fasta     // ppi-splitting's sequences.fasta, keyed by instance id
    family_splits       // channel of [method, split, family-level csv]
    instance_splits     // channel of [method, split, *_instances.csv]

    main:
    instances = INGEST_INSTANCES(domainsplit_db_in, instances_tsv, sequences_fasta)

    family_ordered = orderedSplits(family_splits)
    negatives = INGEST_SAMPLED_NEGATIVES(
        instances.domainsplit_db,
        family_ordered.map { rows -> rows.collect { row -> "${row[0]}:${row[1]}" } },
        family_ordered.map { rows -> rows.collect { row -> row[2] } },
    )

    instance_ordered = orderedSplits(instance_splits)
    membership = INGEST_SPLIT_MEMBERSHIP(
        negatives.domainsplit_db,
        instance_ordered.map { rows -> rows.collect { row -> "${row[0]}:${row[1]}" } },
        instance_ordered.map { rows -> rows.collect { row -> row[2] } },
    )

    ch_versions = Channel.empty().mix(
        INGEST_INSTANCES.out.versions,
        INGEST_SAMPLED_NEGATIVES.out.versions,
        INGEST_SPLIT_MEMBERSHIP.out.versions,
    )

    emit:
    domainsplit_db     = membership.domainsplit_db
    protein_domain_map = instances.protein_domain_map
    unresolved         = membership.unresolved
    versions           = ch_versions
}
