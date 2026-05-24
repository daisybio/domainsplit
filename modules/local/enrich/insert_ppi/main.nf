process INSERT_PPI {
    tag "insert_ppi"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path string_ppi
    path uniprot_id_mapping

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    """
    #!/usr/bin/env python3
    import shutil
    import sqlite3
    import gzip
    import sys
    import pandas as pd
    from tqdm import tqdm
    from typing import Dict

    def load_uniprot_id_mapping(key_name) -> Dict[str, str]:
        mapping = dict()
        with gzip.open("${uniprot_id_mapping}", "rt") as f:
            for line in f:
                uniprot_id, id_type, symbol_name = line.strip().split("\t")
                if id_type.strip() == key_name:
                    mapping[symbol_name.strip()] = uniprot_id.strip()
        return mapping

    shutil.copy("${domainsplit_db_in}", "domainsplit.sqlite3")
    conn = sqlite3.connect("domainsplit.sqlite3")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    print("Inserting PPI", flush=True)

    # Pre-load uniprot_id → protein.id mapping to avoid per-row subqueries
    uniprot_to_pid = dict(
        conn.execute("SELECT uniprot_id, id FROM protein").fetchall()
    )
    print(f"Loaded {len(uniprot_to_pid)} protein ID mappings", flush=True)

    string_id_mapping = load_uniprot_id_mapping("STRING")
    ppi_df = pd.read_csv("${string_ppi}", sep=" ")

    insert_rows = []
    for _, row in tqdm(ppi_df.iterrows(), total=len(ppi_df)):
        uniprot_a = string_id_mapping.get(row["protein1"])
        uniprot_b = string_id_mapping.get(row["protein2"])
        if uniprot_a and uniprot_b:
            pid_a = uniprot_to_pid.get(uniprot_a)
            pid_b = uniprot_to_pid.get(uniprot_b)
            if pid_a is not None and pid_b is not None:
                insert_rows.append((pid_a, pid_b, row["combined_score"]))

    conn.executemany(
        "INSERT INTO protein_protein_interaction(protein_id_a, protein_id_b, score) "
        "VALUES (?, ?, ?)",
        insert_rows,
    )
    conn.commit()
    conn.close()

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
    """
}
