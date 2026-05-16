process SPLIT_DATABASE {
    tag { meta.id }
    label 'process_high'
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path("interactions.csv")
    path "cobinet.sqlite3"

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


    input_db_path = "cobinet.sqlite3"
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
        WITH ddi_to_delete AS (
            SELECT id FROM domain_domain_interaction
            WHERE id NOT IN (SELECT ddi_id FROM database_split)
        )
        DELETE FROM domain_domain_interaction
        WHERE id IN ddi_to_delete;
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

    conn.close()

    import sys as _sys
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {_sys.version.split()[0]}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
    """

}