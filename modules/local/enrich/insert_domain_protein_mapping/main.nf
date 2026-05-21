process INSERT_DOMAIN_PROTEIN_MAPPING {
    tag "insert_domain_protein_mapping"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path protein_domain_map
    path esm_domain_embeddings

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
    import numpy as np
    import h5py
    from tqdm import tqdm

    shutil.copy("${domainsplit_db_in}", "domainsplit.sqlite3")
    conn = sqlite3.connect("domainsplit.sqlite3")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    print("Inserting domain-protein mapping", flush=True)

    def iterate_domain_protein_alignments():
        with gzip.open("${protein_domain_map}", "rt") as protein_domain_map_text, \
            h5py.File("${esm_domain_embeddings}", mode="r") as domain_embeddings:
            pd_map = pd.read_csv(protein_domain_map_text)
            for row in pd_map.itertuples():
                key = f"{row.pfam_id}_{row.uniprot_id}_{row.start_pos}_{row.end_pos}"
                esm3_key = f"{key}/esm3"
                esmc_key = f"{key}/esmc"
                esm3_per_domain = esmc_per_domain = None
                if esm3_key in domain_embeddings:
                    esm3_per_domain = np.array(domain_embeddings[esm3_key]).dumps()
                if esmc_key in domain_embeddings:
                    esmc_per_domain = np.array(domain_embeddings[esmc_key]).dumps()
                yield (row.sequence, row.start_pos, row.end_pos, esm3_per_domain, esmc_per_domain, row.pfam_id, row.uniprot_id)

    conn.executemany(\"""
        INSERT OR IGNORE INTO domain_protein_map(
            domain_id, protein_id, domain_sequence,
            start_pos, end_pos,
            esm3_per_domain, esmc_per_domain
        )
        SELECT domain.id as domain_id, protein.id as protein_id,
            ? as domain_sequence, ? as start_pos, ? as end_pos, ? as esm3_per_domain, ? as esmc_per_domain
        FROM domain, protein
        WHERE
            domain.pfam_id = ? AND
            protein.uniprot_id = ?;
    \""", tqdm(iterate_domain_protein_alignments()))
    conn.commit()
    conn.close()

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
        f.write(f"    numpy: {np.__version__}\\n")
        f.write(f"    h5py: {h5py.__version__}\\n")
    """
}
