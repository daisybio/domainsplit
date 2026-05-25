process RANDOM_DDI_SPLIT {
    tag "random_ddi"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path 'domainsplit.sqlite3'
    val split_fractions  // e.g., [("train", 0.6), ("optimization", 0.2), ("test", 0.2)]

    output:
    path('*.sqlite3'), emit: split_dbs
    val output_split_info, emit: split_info
    path "versions.yml", emit: versions

    script:
    def output_file_fraction_dict = [:]
    output_split_info = []

    split_fractions.each { name, fraction ->
        output_file_fraction_dict["${name}"] = fraction
        output_split_info << ["${name}.sqlite3", name]
    }

    def split_fraction_dict_str = output_file_fraction_dict.collect { k, v -> "'${k}': ${v}" }.join(", ")
    def split_fraction_dict_py = "{" + split_fraction_dict_str + "}"

    """
    #!/usr/bin/env python3

    import os
    os.environ["SQLITE_TMPDIR"] = os.getcwd()

    import sqlite3
    import random
    import shutil

    input_db_path = "domainsplit.sqlite3"
    split_fractions = ${split_fraction_dict_py}

    conn = sqlite3.connect(input_db_path)
    ddi_ids = [row[0] for row in conn.execute("SELECT id FROM domain_domain_interaction")]
    conn.close()

    random.shuffle(ddi_ids)
    total_ddis = len(ddi_ids)

    # Assign DDI IDs to splits
    split_assignments = {}
    items = list(split_fractions.items())
    for i, (split_name, fraction) in enumerate(items):
        is_last = (i == len(items) - 1)
        if is_last:
            split_assignments[split_name] = set(ddi_ids)
            ddi_ids = []
        else:
            count = int(total_ddis * fraction)
            split_assignments[split_name] = set(ddi_ids[:count])
            ddi_ids = ddi_ids[count:]

    # Create per-split databases
    for split_name, keep_ids in split_assignments.items():
        output_path = f"{split_name}.sqlite3"
        print(f"Creating {output_path} with {len(keep_ids)} DDIs...")
        shutil.copyfile(input_db_path, output_path)

        conn = sqlite3.connect(output_path)
        conn.executescript('''
            PRAGMA foreign_keys=ON;
            PRAGMA journal_mode=OFF;
            PRAGMA synchronous=OFF;
        ''')

        # Create temp table with IDs to keep
        conn.execute("CREATE TEMP TABLE keep_ids (id INTEGER PRIMARY KEY)")
        conn.executemany("INSERT INTO keep_ids VALUES (?)", [(i,) for i in keep_ids])

        conn.execute('''
            DELETE FROM domain_domain_interaction
            WHERE id NOT IN (SELECT id FROM keep_ids)
        ''')

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

        conn.execute("DROP TABLE keep_ids")

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
        print(f"  {output_path}: done")

    import sys as _sys
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {_sys.version.split()[0]}\\n")
    """
}
