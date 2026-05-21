process INSERT_PROTEINS_WITH_EMBEDDINGS {
    tag "insert_proteins_with_embeddings"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path uniprot_database
    path protein_domain_map
    path prott5_embeddings
    path esm_protein_embeddings

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
    from Bio import SeqIO
    from tqdm import tqdm

    shutil.copy("${domainsplit_db_in}", "domainsplit.sqlite3")
    conn = sqlite3.connect("domainsplit.sqlite3")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    print("Inserting required uniprot records and embeddings", flush=True)
    with gzip.open("${uniprot_database}","rt") as uniprot_text, \
            gzip.open("${protein_domain_map}","rt") as protein_domain_map_text, \
            h5py.File("${esm_protein_embeddings}", mode="r") as esm_protein_embeddings, \
            h5py.File("${prott5_embeddings}", mode="r") as prott5_embeddings_file:

        pd_map = pd.read_csv(protein_domain_map_text)
        required_uniprot_ids = set(pd_map["uniprot_id"].unique())

        uniprot_records = SeqIO.parse(uniprot_text, "fasta")
        uniprot_records = map(lambda record: (record.id.split("|")[1], str(record.seq)), uniprot_records)
        uniprot_records = filter(lambda record: record[0] in required_uniprot_ids, uniprot_records)

        get_prott5_embedding = lambda seq_id: np.array(prott5_embeddings_file[seq_id]).dumps() \
            if seq_id in prott5_embeddings_file else None
        uniprot_records = map(lambda record: record + (get_prott5_embedding(record[0]),), uniprot_records)

        def get_esm_embeddings(record):
            seq_id = record[0]
            esmc_key = f"{seq_id}/esmc"
            esm3_key = f"{seq_id}/esm3"
            esm3_embedding = np.array(esm_protein_embeddings[esm3_key]).dumps() if esm3_key in esm_protein_embeddings else None
            esmc_embedding = np.array(esm_protein_embeddings[esmc_key]).dumps() if esmc_key in esm_protein_embeddings else None
            return (esm3_embedding, esmc_embedding)

        uniprot_records = map(lambda record: record + get_esm_embeddings(record), uniprot_records)
        uniprot_records = tqdm(uniprot_records)

        conn.executemany('''INSERT INTO protein (
            uniprot_id, sequence,
            prott5_per_residue, esm3_per_residue, esmc_per_residue
        )
        VALUES (?, ?, ?, ?, ?);''', uniprot_records)
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
