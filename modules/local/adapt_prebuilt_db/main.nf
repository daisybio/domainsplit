process ADAPT_PREBUILT_DB {
    tag "adapt_prebuilt_db"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker.io/konstantinpelz/domainsplit-general:1.0.0"

    input:
    path db

    output:
    path "domainsplit_adapted.sqlite3", emit: domainsplit_db
    path "versions.yml",                emit: versions

    script:
    """
    #!/usr/bin/env python3
    import shutil
    import sqlite3
    import sys

    src_db = "${db}"
    dest_db = "domainsplit_adapted.sqlite3"


    shutil.copyfile(src_db, dest_db)

    con = sqlite3.connect(dest_db)

    def tables():
        return {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}

    def columns(table):
        return [r[1] for r in con.execute(f"PRAGMA table_info({table})")]

    # Tables the downstream steps need and this script does not create.
    required = ["domain", "domain_domain_interaction", "protein",
                "domain_protein_map", "ddi_split_membership"]
    missing_tables = [t for t in required if t not in tables()]
    if missing_tables:
        sys.exit(f"ERROR: --prebuilt_db lacks required tables: {missing_tables}")

    # Adding missing columns on ddi_split_membership
    wanted = [
        ("z_score", "REAL"),
        ("is_mock", "INTEGER DEFAULT 0"),
    ]
    present = set(columns("ddi_split_membership"))
    for name, ddl in wanted:
        if name in present:
            print(f"ddi_split_membership.{name}: already present")
        else:
            con.execute(f"ALTER TABLE ddi_split_membership ADD COLUMN {name} {ddl}")
            print(f"ddi_split_membership.{name}: added")

    # Adding missing table ddi_interaction_confirmed
    con.execute('''
        CREATE TABLE IF NOT EXISTS ddi_interaction_confirmed (
            ddi_id REFERENCES domain_domain_interaction ON DELETE CASCADE,
            method,
            majority_confirmed INTEGER,
            mean_confirmed INTEGER,
            UNIQUE(ddi_id, method)
        )
    ''')

    # Adding missing indexes
    con.executescript('''
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
        CREATE INDEX IF NOT EXISTS idx_ddi_split_membership_triple
        ON ddi_split_membership (ddi_id, instance_id_a, instance_id_b);
    ''')
    con.commit()

    # Print checks
    cols = columns("ddi_split_membership")
    print("ddi_split_membership columns:", cols)
    for t in ["ddi_split_membership", "ddi_interaction_confirmed", "domain_protein_map"]:
        n = con.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0]
        print(f"{t}: {n} rows")
    for m, s, n in con.execute(
        "SELECT method, split, COUNT(*) FROM ddi_split_membership GROUP BY method, split"
    ):
        print(f"  membership {m}/{s}: {n}")
    con.close()

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\\n")
    """
}