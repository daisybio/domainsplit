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

    def ppi_iterator():
        string_id_mapping = load_uniprot_id_mapping("STRING")
        ppi_df = pd.read_csv("${string_ppi}", sep=" ")

        for _, row in ppi_df.iterrows():
            uniprot_id_a = string_id_mapping.get(row["protein1"])
            uniprot_id_b = string_id_mapping.get(row["protein2"])
            if uniprot_id_a and uniprot_id_b:
                yield row["combined_score"], uniprot_id_a, uniprot_id_b

    conn.executemany(\"""
        INSERT INTO protein_protein_interaction(protein_id_a, protein_id_b, score)
        SELECT protein_a.id, protein_b.id, ? as score
        FROM protein AS protein_a, protein AS protein_b
        WHERE
            protein_a.uniprot_id = ? AND
            protein_b.uniprot_id = ?;
    \""", tqdm(ppi_iterator()))
    conn.commit()
    conn.close()

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
    """
}
