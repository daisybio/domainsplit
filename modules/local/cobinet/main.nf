process COBINET {
    tag "cobinet"
    label 'process_high'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path ddi_3did
    path ddi_negatome
    path pfam2go
    path uniprot_database
    path protein_domain_map
    path prott5_embeddings
    path uniprot_go_terms
    path string_ppi
    path uniprot_id_mapping
    path esm_protein_embeddings
    path esm_domain_embeddings

    output:
    path "cobinet.sqlite3", emit: cobinet_db
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3

    import sqlite3
    import gzip
    import itertools
    import pandas as pd
    import numpy as np
    import h5py
    import re
    from Bio import SeqIO
    from Bio import AlignIO
    from tqdm import tqdm
    from typing import Set, Dict
    from pathlib import Path

    def load_pfam2go_dataframe(pfam2go_path: str) -> pd.DataFrame:
        pfam2go_text = Path(pfam2go_path).read_text()
        matches = re.finditer(
            r"^Pfam:(?P<Pfam_accession>PF\\d*)\\s*(?P<Pfam_name>.*?)\\s*>\\s*GO:(?P<GO_name>.*?)\\s*;\\s*(?P<GO_accession>.*)\$",
            pfam2go_text,
            flags=re.MULTILINE
        )
        return pd.DataFrame.from_records(match.groupdict() for match in matches)

    def load_uniprot_id_mapping(key_name) -> Dict[str, str]:
        mapping = dict()
        with gzip.open("${uniprot_id_mapping}","rt") as f:
            for line in f:
                uniprot_id, id_type, symbol_name = line.strip().split("\t")
                if id_type.strip() == key_name:
                    mapping[symbol_name.strip()] = uniprot_id.strip()
        return mapping

    def get_negatome_pfam_ids() -> Set[str]:
        pfam_ids = set()
        with open("${ddi_negatome}") as f:
            for line in f:
                id_1, id_2 = line.strip().split()
                pfam_ids.add(id_1)
                pfam_ids.add(id_2)
        return pfam_ids

    conn_3did = sqlite3.connect("${ddi_3did}")
    conn_cobinet = sqlite3.connect("cobinet.sqlite3")
    print("Creating database schema")
    conn_cobinet.executescript(\"""
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
        \""")

    # SELECT name and Pfam id from 3did and move to cobinet database
    cursor = conn_3did.execute(\"""
        SELECT Name, Pfam_id, profile_length FROM Domain, domain_length
        WHERE domain_length.domain = Domain.Name
        \""")

    print("Inserting domain information")
    insert_iterator_3did = ((name, pfam_id.split(".")[0]) for (name, pfam_id, length) in cursor)
    insert_iterator_negatome = ((None, pfam_id) for pfam_id in get_negatome_pfam_ids())

    insert_iterator = tqdm(itertools.chain(insert_iterator_3did,insert_iterator_negatome))

    # need to insert or ignore, because the same pfam id can occur in both 3did and negatome. This eliminates
    # duplicates that would violate the unique constraint
    conn_cobinet.executemany(
        "INSERT OR IGNORE INTO domain(name, pfam_id) VALUES (?, ?);",
        insert_iterator)
    cursor.close()
    conn_cobinet.commit()

    print("Inserting domain-domain interactions")
    # SELECT the domain domain interactions and tranfer over to cobinet
    cursor = conn_3did.execute(\"""
                SELECT d1.Pfam_id, d2.Pfam_id
                FROM DDI1, Domain AS d1, Domain AS d2
                WHERE
                    DDI1.domain1 = d1.Name AND
                    DDI1.domain2 = d2.Name;
                \""")

    pfam_id_iterator = tqdm((id_1.split(".")[0], id_2.split(".")[0]) for (id_1, id_2) in cursor)
    conn_cobinet.executemany(\"""
                    INSERT INTO domain_domain_interaction(domain_id_a, domain_id_b, negative)
                    SELECT d1.id AS domain_id_a, d2.id AS domain_id_b, FALSE as negative
                    FROM domain AS d1, domain AS d2
                    WHERE d1.pfam_id = ? AND d2.pfam_id = ?
                    ;
                \""",pfam_id_iterator)
    cursor.close()
    conn_cobinet.commit()

    print("Inserting negative domain-domain intereactions (negatome)")
    with open("${ddi_negatome}") as f:
        pfam_id_iterator = ((line.split() for line in f))
        try:
            # TODO investigate why negatome contains some erroneous entries
            conn_cobinet.executemany(\"""
                            INSERT OR IGNORE INTO domain_domain_interaction(domain_id_a, domain_id_b, negative)
                            SELECT d1.id AS domain_id_a, d2.id AS domain_id_b, TRUE as negative
                            FROM domain AS d1, domain AS d2
                            WHERE d1.pfam_id = ? AND d2.pfam_id = ?
                            ;
                        \""",pfam_id_iterator)
        except Exception as e:
            print(next(pfam_id_iterator))
            raise e
    conn_cobinet.commit()

    print("Inserting domain GO information")
    pfam2go_dataframe = load_pfam2go_dataframe("${pfam2go}")
    insert_iterator = ((row["GO_accession"], row["Pfam_accession"]) for i, row in pfam2go_dataframe.iterrows())
    conn_cobinet.executemany(\"""
        INSERT INTO domain_go_terms(domain_id, go_accession)
        SELECT domain.id as domain_id, ? AS go_accession
        FROM domain
        WHERE domain.pfam_id = ?;
    \""",insert_iterator)
    conn_cobinet.commit()

    print("Inserting required uniprot records and embeddings")
    # filter the uniprot fasta and put into cobinet database
    with gzip.open("${uniprot_database}","rt") as uniprot_text, \
            gzip.open("${protein_domain_map}","rt") as protein_domain_map_text, \
            h5py.File("${esm_protein_embeddings}",mode="r") as esm_protein_embeddings, \
            h5py.File("${prott5_embeddings}",mode="r") as prott5_embeddings_file:

        pd_map = pd.read_csv(protein_domain_map_text)
        required_uniprot_ids = set(pd_map["uniprot_id"].unique())

        uniprot_records = SeqIO.parse(uniprot_text,"fasta")
        # convert uniprot records to (id, seq) tuples
        uniprot_records = map(lambda record: (record.id.split("|")[1], str(record.seq)),uniprot_records)

        # filter for required records
        uniprot_records = filter(lambda record: record[0] in required_uniprot_ids,uniprot_records)

        # add embeddings to tuple
        get_prott5_embedding = lambda seq_id: np.array(prott5_embeddings_file[seq_id]).dumps() \
            if seq_id in prott5_embeddings_file else None
        uniprot_records = map(lambda record: record + (get_prott5_embedding(record[0]),),uniprot_records)

        def get_esm_embeddings(record):
            seq_id = record[0]
            esmc_key = f"{seq_id}/esmc"
            esm3_key = f"{seq_id}/esm3"
            esm3_embedding = np.array(esm_protein_embeddings[esm3_key]).dumps() if esm3_key in esm_protein_embeddings else None
            esmc_embedding = np.array(esm_protein_embeddings[esmc_key]).dumps() if esmc_key in esm_protein_embeddings else None
            return (esm3_embedding, esmc_embedding)

        uniprot_records = map(lambda record: record + get_esm_embeddings(record),uniprot_records)

        # filter out records without embeddings
        #uniprot_records = filter(lambda record: record[2] is not None,uniprot_records)

        # tqdm statistics
        uniprot_records = tqdm(uniprot_records)

        conn_cobinet.executemany('''INSERT INTO protein (
            uniprot_id, sequence,
            prott5_per_residue, esm3_per_residue, esmc_per_residue
        )
        VALUES (?, ?, ?, ?, ?);''', uniprot_records)
    conn_cobinet.commit()

    print("Inserting protein GO information")
    # read the uniprot go terms and put them into the database
    with gzip.open("${uniprot_go_terms}","rt") as go_terms_text:
        go_terms_df = pd.read_csv(go_terms_text,sep="\t")
        # Proteins without GO annotation are NaN; skip them before splitting.
        insert_iterator = (((go_term, row["Entry"]) for go_term in row["Gene Ontology IDs"].split("; "))
                           for _, row in go_terms_df.iterrows()
                           if isinstance(row["Gene Ontology IDs"], str))
        insert_iterator = itertools.chain.from_iterable(insert_iterator)
        insert_iterator = tqdm(insert_iterator)

        conn_cobinet.executemany(\"""
            INSERT INTO protein_go_terms(protein_id, go_accession)
            SELECT protein.id as protein_id, ? AS go_accession
            FROM protein
            WHERE protein.uniprot_id = ?;
        \""",insert_iterator)

    print("Inserting PPI")


    def ppi_iterator():
        string_id_mapping = load_uniprot_id_mapping("STRING")
        ppi_df = pd.read_csv("${string_ppi}",sep=" ")

        for _, row in ppi_df.iterrows():
            uniprot_id_a = string_id_mapping.get(row["protein1"])
            uniprot_id_b = string_id_mapping.get(row["protein2"])
            if uniprot_id_a and uniprot_id_b:
                yield row["combined_score"], uniprot_id_a, uniprot_id_b


    conn_cobinet.executemany(\"""
        INSERT INTO protein_protein_interaction(protein_id_a, protein_id_b, score)
        SELECT protein_a.id, protein_b.id, ? as score
        FROM protein AS protein_a, protein AS protein_b
        WHERE
            protein_a.uniprot_id = ? AND
            protein_b.uniprot_id = ?;
    \""",tqdm(ppi_iterator()))
    conn_cobinet.commit()

    # read alignments and put domain-protein mappings into cobinet
    print("Inserting domain-protein mapping")

    def iterate_domain_protein_alignments():
        with gzip.open("${protein_domain_map}","rt") as protein_domain_map_text, \
            h5py.File("${esm_domain_embeddings}",mode="r") as domain_embeddings:
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

    conn_cobinet.executemany(\"""
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
    \""",tqdm(iterate_domain_protein_alignments()))

    conn_cobinet.commit()
    conn_cobinet.close()

    import sys as _sys
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {_sys.version.split()[0]}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
        f.write(f"    numpy: {np.__version__}\\n")
        f.write(f"    h5py: {h5py.__version__}\\n")
    """
}
