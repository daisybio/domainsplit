/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBSET_DDIS_BY_SOURCE -- build a single split database keeping only DDIs
    whose `source` is in the requested set, then prune orphan domains/proteins.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Used for the External-Validation test set, which is placed "as is" (no
    leakage-aware partitioning) from the held-out sources.
----------------------------------------------------------------------------*/

process SUBSET_DDIS_BY_SOURCE {
    tag "subset_${split_name}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path 'domainsplit.sqlite3'
    val  source_filter   // list of DDI source strings to keep
    val  split_name      // output split name, e.g. 'test'

    output:
    path('*.sqlite3'), emit: split_dbs
    val output_split_info, emit: split_info
    path "versions.yml", emit: versions

    script:
    output_split_info = [["${split_name}.sqlite3", split_name]]
    def src_list = source_filter.collect { "'${it}'" }.join(", ")

    """
    #!/usr/bin/env python3
    import os
    os.environ["SQLITE_TMPDIR"] = os.getcwd()

    import sqlite3
    import shutil
    import sys

    input_db_path = "domainsplit.sqlite3"
    output_path = "${split_name}.sqlite3"
    sources = (${src_list},)

    shutil.copyfile(input_db_path, output_path)

    conn = sqlite3.connect(output_path)
    conn.executescript('''
        PRAGMA foreign_keys=ON;
        PRAGMA journal_mode=OFF;
        PRAGMA synchronous=OFF;
    ''')

    placeholders = ",".join("?" for _ in sources)
    n_keep = conn.execute(
        f"SELECT COUNT(*) FROM domain_domain_interaction WHERE source IN ({placeholders})",
        sources,
    ).fetchone()[0]
    print(f"Keeping {n_keep} DDIs with source in {sources}", flush=True)

    conn.execute(
        f"DELETE FROM domain_domain_interaction WHERE source NOT IN ({placeholders})",
        sources,
    )

    conn.execute('''
        DELETE FROM domain WHERE id IN (
            SELECT d.id FROM domain d
            LEFT JOIN domain_domain_interaction ddi
                ON ddi.domain_id_a = d.id OR ddi.domain_id_b = d.id
            LEFT JOIN domain_protein_map dpm
                ON dpm.domain_id = d.id
            WHERE ddi.id IS NULL OR dpm.domain_id IS NULL
        )
    ''')

    conn.execute('''
        DELETE FROM protein WHERE id IN (
            SELECT p.id FROM protein p
            LEFT JOIN domain_protein_map dpm
                ON dpm.protein_id = p.id
            WHERE dpm.domain_id IS NULL
        )
    ''')

    conn.executescript('''
        VACUUM;

        CREATE INDEX IF NOT EXISTS idx_ddi_domain_a ON domain_domain_interaction(domain_id_a);
        CREATE INDEX IF NOT EXISTS idx_ddi_domain_b ON domain_domain_interaction(domain_id_b);
        CREATE INDEX IF NOT EXISTS idx_dpm_domain ON domain_protein_map(domain_id);
        CREATE INDEX IF NOT EXISTS idx_dpm_protein ON domain_protein_map(protein_id);
        CREATE INDEX IF NOT EXISTS idx_ppi_protein_a ON protein_protein_interaction(protein_id_a);
        CREATE INDEX IF NOT EXISTS idx_ppi_protein_b ON protein_protein_interaction(protein_id_b);
        CREATE INDEX IF NOT EXISTS idx_pgo_protein ON protein_go_terms(protein_id);
    ''')

    conn.close()
    print(f"  {output_path}: done", flush=True)

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\\n")
    """
}
