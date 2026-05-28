process INIT_DOMAINSPLIT_DB {
    tag "init_domainsplit_db"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

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
        CREATE TABLE domain_domain_interaction (
            id INTEGER PRIMARY KEY,
            domain_id_a, domain_id_b, negative,
            source VARCHAR(255),
            FOREIGN KEY(domain_id_a) REFERENCES domain ON DELETE CASCADE,
            FOREIGN KEY(domain_id_b) REFERENCES domain ON DELETE CASCADE,
            UNIQUE(domain_id_a, domain_id_b)
        );

        CREATE TABLE protein (
            id INTEGER PRIMARY KEY,
            uniprot_id,
            sequence,
            prott5_per_residue,
            esm3_per_residue,
            esmc_per_residue,
            UNIQUE(uniprot_id)
        );
        CREATE TABLE protein_go_terms(
            protein_id REFERENCES protein ON DELETE CASCADE,
            go_accession
        );
        CREATE TABLE protein_protein_interaction (
            protein_id_a REFERENCES protein ON DELETE CASCADE,
            protein_id_b REFERENCES protein ON DELETE CASCADE,
            score,
            UNIQUE(protein_id_a, protein_id_b)
        );

        CREATE TABLE domain_protein_map (
            domain_id REFERENCES domain ON DELETE CASCADE,
            protein_id REFERENCES protein ON DELETE CASCADE,
            domain_sequence, start_pos, end_pos,
            esm3_per_domain, esmc_per_domain,
            UNIQUE(domain_id, protein_id)
        );

        CREATE INDEX IF NOT EXISTS idx_domain_domain_interaction_domain_id_a
        ON domain_domain_interaction (domain_id_a);
        CREATE INDEX IF NOT EXISTS idx_domain_domain_interaction_domain_id_b
        ON domain_domain_interaction (domain_id_b);
        CREATE INDEX IF NOT EXISTS idx_domain_protein_map_domain_id
        ON domain_protein_map (domain_id);
        CREATE INDEX IF NOT EXISTS idx_domain_protein_map_protein_id
        ON domain_protein_map (protein_id);
        CREATE INDEX IF NOT EXISTS idx_protein_protein_interaction_protein_id_a
        ON protein_protein_interaction (protein_id_a);
        CREATE INDEX IF NOT EXISTS idx_protein_protein_interaction_protein_id_b
        ON protein_protein_interaction (protein_id_b);
    ''')
    con.commit()
    con.close()

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\\n")
    """
}
