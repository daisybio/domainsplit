process SPLIT_DATABASE {
    tag { meta.id }
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    tuple val(meta), path("interactions.csv")
    path "domainsplit.sqlite3"

    output:
    tuple val(meta), path(output_db_path), emit: split_db
    path "versions.yml", emit: versions

    script:
    output_db_path = "${meta.id}.sqlite3"

    """
    #!/usr/bin/env python3

    import os

    # set SQLITE_TMPDIR to current working directory to avoid issues with temp files on some filesystems
    os.environ["SQLITE_TMPDIR"] = os.getcwd()

    import sqlite3
    import shutil
    import pandas as pd


    input_db_path = "domainsplit.sqlite3"
    output_db_path = "${output_db_path}"
    interaction_ids_path = "interactions.csv"

    # Copy the database to the output path
    shutil.copyfile(input_db_path, output_db_path)

    # Connect to the copied database
    conn = sqlite3.connect(output_db_path)

    # Read interaction IDs from the provided text file
    # and input them into a temporary table for faster lookup
    print("Creating temporary table for split interaction IDs...")
    conn.executescript('''
        PRAGMA foreign_keys=ON;
        PRAGMA journal_mode=OFF;
        PRAGMA synchronous=OFF;

        CREATE TABLE database_split (
            ddi_id UNIQUE REFERENCES domain_domain_interaction ON DELETE CASCADE,
            is_evaluation_relevant BOOL DEFAULT 1,
            override_negative BOOL DEFAULT NULL
        );

        ALTER TABLE domain_domain_interaction ADD COLUMN is_evaluation_relevant;
    ''')

    print("Inserting interaction IDs into temporary table...")
    interactions_df = pd.read_csv(interaction_ids_path)
    interactions_df.to_sql('database_split', conn, if_exists='append', index=False)

    # Delete entries from domain_domain_interaction not in database_split
    print("Deleting non-split entries from domain_domain_interaction...")
    conn.executescript('''
        DELETE FROM domain_domain_interaction
        WHERE NOT EXISTS (
            SELECT 1 FROM database_split
            WHERE database_split.ddi_id = domain_domain_interaction.id
        );
    ''')

    print("Updating remaining domain_domain_interaction entries")
    conn.executescript('''
        UPDATE domain_domain_interaction
        SET is_evaluation_relevant = split.is_evaluation_relevant,
            negative = COALESCE(split.override_negative, negative)
        FROM database_split AS split
        WHERE split.ddi_id = id;
    ''')

    print("Removing orphaned domains...")
    # orphaned means not referenced in domain_domain_interaction or
    # missing from domain_protein_map
    conn.execute('''
        WITH orphaned_domains AS (
            SELECT domain.id FROM domain
            LEFT OUTER JOIN domain_domain_interaction
            ON domain_domain_interaction.domain_id_a = domain.id OR
                domain_domain_interaction.domain_id_b = domain.id
            LEFT OUTER JOIN domain_protein_map
            ON domain_protein_map.domain_id = domain.id
            WHERE domain_domain_interaction.id IS NULL
                OR domain_protein_map.domain_id IS NULL
        )
        DELETE FROM domain WHERE id IN orphaned_domains;
    ''')

    print("Removing orphaned proteins...")
    conn.execute('''
        WITH orphaned_proteins AS (
            SELECT protein.id FROM protein
            LEFT OUTER JOIN domain_protein_map
            ON domain_protein_map.protein_id = protein.id
            WHERE domain_protein_map.domain_id IS NULL
        )
        DELETE FROM protein WHERE id IN orphaned_proteins;
    ''')

    print("Vacuuming the database to optimize size...")
    conn.executescript('''
        DROP TABLE database_split;
        VACUUM;
    ''')

    print("Creating indexes for downstream benchmark queries...")
    conn.executescript('''
        CREATE INDEX IF NOT EXISTS idx_ddi_domain_a ON domain_domain_interaction(domain_id_a);
        CREATE INDEX IF NOT EXISTS idx_ddi_domain_b ON domain_domain_interaction(domain_id_b);
        CREATE INDEX IF NOT EXISTS idx_ddi_eval_relevant ON domain_domain_interaction(is_evaluation_relevant);
        CREATE INDEX IF NOT EXISTS idx_dpm_domain ON domain_protein_map(domain_id);
        CREATE INDEX IF NOT EXISTS idx_dpm_protein ON domain_protein_map(protein_id);
        CREATE INDEX IF NOT EXISTS idx_ppi_protein_a ON protein_protein_interaction(protein_id_a);
        CREATE INDEX IF NOT EXISTS idx_ppi_protein_b ON protein_protein_interaction(protein_id_b);
        CREATE INDEX IF NOT EXISTS idx_pgo_protein ON protein_go_terms(protein_id);
    ''')

    conn.close()

    import sys as _sys
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {_sys.version.split()[0]}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
    """

}
