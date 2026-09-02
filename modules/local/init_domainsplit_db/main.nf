process INIT_DOMAINSPLIT_DB {
    tag "init_domainsplit_db"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    """
    #!/usr/bin/env python3
    import sqlite3
    import sys

    con = sqlite3.connect("domainsplit.sqlite3")
    con.executescript('''
        PRAGMA foreign_keys=ON;
        PRAGMA journal_mode=OFF;
        PRAGMA synchronous=OFF;

        CREATE TABLE domain (id INTEGER PRIMARY KEY, pfam_id, name, UNIQUE(pfam_id));
        CREATE TABLE domain_go_terms(
            domain_id REFERENCES domain ON DELETE CASCADE,
            go_accession
        );
        -- One row per domain pair. `source` is the comma-joined list of every
        -- source that contributed the pair (e.g. "PPIDM,PPIDM_Gold"); see
        -- bin/ddi_db_utils.py for the merge/drop rules that maintain it.
        CREATE TABLE domain_domain_interaction (
            id INTEGER PRIMARY KEY,
            domain_id_a, domain_id_b, negative,
            source VARCHAR(255),
            FOREIGN KEY(domain_id_a) REFERENCES domain ON DELETE CASCADE,
            FOREIGN KEY(domain_id_b) REFERENCES domain ON DELETE CASCADE,
            UNIQUE(domain_id_a, domain_id_b)
        );

        -- `taxon_id` and `reviewed` both come from instances.tsv via
        -- INGEST_INSTANCES (`reviewed` out of its `source_db` column). They are
        -- one value each on a `human_reviewed` run, but without them a
        -- multi-species / multi-review-status database cannot be stratified, and
        -- there is no way to tell which proteins went unenriched because their
        -- species had no STRING file or because UniProt cross-references no STRING
        -- id for unreviewed entries.
        CREATE TABLE protein (
            id INTEGER PRIMARY KEY,
            uniprot_id,
            sequence,
            taxon_id,
            reviewed TEXT,
            UNIQUE(uniprot_id)
        );
        CREATE TABLE protein_go_terms(
            protein_id REFERENCES protein ON DELETE CASCADE,
            go_accession
        );
        -- `score` is typed REAL for the same reason start_pos/end_pos below are
        -- typed INTEGER: declared with no type the column takes BLOB (none)
        -- affinity, so SQLite stores whatever class the writer binds. INSERT_PPI
        -- used to bind STRING's combined_score as the raw string it split out of
        -- the links file, which made every reader see TEXT '400' -- and a
        -- consumer comparing it to a numeric confidence cutoff ('score >= 400')
        -- raised instead of filtering. INSERT_PPI parses the score now; the
        -- affinity is the second guard, so a writer binding '400' is coerced to
        -- 400.0 and no reader has to guess.
        CREATE TABLE protein_protein_interaction (
            protein_id_a REFERENCES protein ON DELETE CASCADE,
            protein_id_b REFERENCES protein ON DELETE CASCADE,
            score REAL,
            UNIQUE(protein_id_a, protein_id_b)
        );

        -- One row per domain *instance*, matching the embedding H5 key contract
        -- ({pfam_id}_{uniprot_id}_{start}_{end}): a protein carrying two copies
        -- of the same family gets two rows. `instance_id` is ppi-splitting's
        -- own instance identifier (NULL for rows not sourced from it).
        --
        -- start_pos/end_pos are typed INTEGER on purpose. Declared with no type
        -- they take BLOB (none) affinity, SQLite converts nothing on insert, and
        -- TEXT '10' is not INTEGER 10 in the UNIQUE key below -- which is how two
        -- writers that disagreed on the storage class ended up producing two rows
        -- per instance instead of upserting one. With the affinity declared, a
        -- writer binding '10' is coerced to 10 and cannot desync again.
        CREATE TABLE domain_protein_map (
            domain_id REFERENCES domain ON DELETE CASCADE,
            protein_id REFERENCES protein ON DELETE CASCADE,
            domain_sequence, start_pos INTEGER, end_pos INTEGER,
            instance_id, clan, taxon_id,
            UNIQUE(domain_id, protein_id, start_pos, end_pos)
        );

        -- Which DDIs, and which concrete domain-instance pairs, belong to which
        -- (method, split). Populated by INGEST_SPLIT_MEMBERSHIP from
        -- ppi-splitting's instance-level CSVs and by BUILD_EXTERNAL_TEST.
        -- Exists only to drive SUBSET_SPLIT_DB, which is a pure SQL filter over
        -- it. `instance_id_*` reference domain_protein_map.instance_id rather
        -- than protein.id: a protein carrying two copies of one family gives two
        -- instances, and the split is defined on instances, not proteins.
        CREATE TABLE ddi_split_membership (
            ddi_id REFERENCES domain_domain_interaction ON DELETE CASCADE,
            method, split,
            instance_id_a, instance_id_b,
            UNIQUE(ddi_id, method, split, instance_id_a, instance_id_b)
        );

        CREATE INDEX IF NOT EXISTS idx_domain_domain_interaction_domain_id_a
        ON domain_domain_interaction (domain_id_a);
        CREATE INDEX IF NOT EXISTS idx_domain_domain_interaction_domain_id_b
        ON domain_domain_interaction (domain_id_b);
        CREATE INDEX IF NOT EXISTS idx_domain_protein_map_domain_id
        ON domain_protein_map (domain_id);
        CREATE INDEX IF NOT EXISTS idx_domain_protein_map_protein_id
        ON domain_protein_map (protein_id);
        CREATE UNIQUE INDEX IF NOT EXISTS idx_domain_protein_map_instance_id
        ON domain_protein_map (instance_id) WHERE instance_id IS NOT NULL;
        CREATE INDEX IF NOT EXISTS idx_protein_protein_interaction_protein_id_a
        ON protein_protein_interaction (protein_id_a);
        CREATE INDEX IF NOT EXISTS idx_protein_protein_interaction_protein_id_b
        ON protein_protein_interaction (protein_id_b);
        CREATE INDEX IF NOT EXISTS idx_ddi_split_membership_method_split
        ON ddi_split_membership (method, split);
        CREATE INDEX IF NOT EXISTS idx_ddi_split_membership_ddi_id
        ON ddi_split_membership (ddi_id);
    ''')
    con.commit()
    con.close()

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\\n")
    """
}
