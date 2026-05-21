process INSERT_PROTEIN_GO_TERMS {
    tag "insert_protein_go_terms"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path uniprot_go_terms

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    """
    #!/usr/bin/env python3
    import shutil
    import sqlite3
    import gzip
    import itertools
    import sys
    import pandas as pd
    from tqdm import tqdm

    shutil.copy("${domainsplit_db_in}", "domainsplit.sqlite3")
    conn = sqlite3.connect("domainsplit.sqlite3")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    print("Inserting protein GO information", flush=True)
    with gzip.open("${uniprot_go_terms}","rt") as go_terms_text:
        go_terms_df = pd.read_csv(go_terms_text, sep="\t")
        # Proteins without GO annotation are NaN; skip them before splitting.
        insert_iterator = (((go_term, row["Entry"]) for go_term in row["Gene Ontology IDs"].split("; "))
                           for _, row in go_terms_df.iterrows()
                           if isinstance(row["Gene Ontology IDs"], str))
        insert_iterator = itertools.chain.from_iterable(insert_iterator)
        insert_iterator = tqdm(insert_iterator)

        conn.executemany(\"""
            INSERT INTO protein_go_terms(protein_id, go_accession)
            SELECT protein.id as protein_id, ? AS go_accession
            FROM protein
            WHERE protein.uniprot_id = ?;
        \""", insert_iterator)
    conn.commit()
    conn.close()

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
    """
}
