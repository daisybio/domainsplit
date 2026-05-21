process INSERT_DOMAIN_GO_TERMS {
    tag "insert_domain_go_terms"
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path domainsplit_db_in, stageAs: 'input.domainsplit.sqlite3'
    path pfam2go

    output:
    path "domainsplit.sqlite3", emit: domainsplit_db
    path "versions.yml",        emit: versions

    script:
    """
    #!/usr/bin/env python3
    import shutil
    import sqlite3
    import re
    import sys
    import pandas as pd
    from pathlib import Path

    def load_pfam2go_dataframe(pfam2go_path: str) -> pd.DataFrame:
        pfam2go_text = Path(pfam2go_path).read_text()
        matches = re.finditer(
            r"^Pfam:(?P<Pfam_accession>PF\\d*)\\s*(?P<Pfam_name>.*?)\\s*>\\s*GO:(?P<GO_name>.*?)\\s*;\\s*(?P<GO_accession>.*)\$",
            pfam2go_text,
            flags=re.MULTILINE
        )
        return pd.DataFrame.from_records(match.groupdict() for match in matches)

    shutil.copy("${domainsplit_db_in}", "domainsplit.sqlite3")
    conn = sqlite3.connect("domainsplit.sqlite3")
    conn.execute("PRAGMA foreign_keys=ON")
    conn.execute("PRAGMA journal_mode=OFF")
    conn.execute("PRAGMA synchronous=OFF")

    print("Inserting domain GO information", flush=True)
    pfam2go_df = load_pfam2go_dataframe("${pfam2go}")
    insert_iterator = ((row["GO_accession"], row["Pfam_accession"]) for _, row in pfam2go_df.iterrows())
    conn.executemany(\"""
        INSERT INTO domain_go_terms(domain_id, go_accession)
        SELECT domain.id as domain_id, ? AS go_accession
        FROM domain
        WHERE domain.pfam_id = ?;
    \""", insert_iterator)
    conn.commit()
    conn.close()

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    pandas: {pd.__version__}\\n")
    """
}
