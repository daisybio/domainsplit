process RANDOM_DDI_SPLIT {
    tag "random_ddi"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"

    input:
    path 'cobinet.sqlite3'
    val split_fractions  // e.g., [("train", 0.8), ("test", 0.2)]

    output:
    path('*.csv'), emit: split_ddi_id_files
    val output_file_splits, emit: split_fractions
    path "versions.yml", emit: versions

    script:
    def output_file_fraction_dict = [:]
    output_file_splits = []

    split_fractions.each { name, fraction ->
        output_file_fraction_dict["${name}.csv"] = fraction
        output_file_splits << ["${name}.csv", name]
    }

    def split_fraction_dict_str = output_file_fraction_dict.collect { k, v -> "'${k}': ${v}" }.join(", ")
    def split_fraction_dict_py = "{" + split_fraction_dict_str + "}"

    """
    #!/usr/bin/env python3

    import sqlite3
    import random
    import os

    input_db_path = "cobinet.sqlite3"
    output_file_fraction_dict = ${split_fraction_dict_py}

    # Connect to the database
    conn = sqlite3.connect(input_db_path)
    # Fetch all domain-domain interaction IDs
    ddi_ids = conn.execute("SELECT id FROM domain_domain_interaction;")
    ddi_ids = [row[0] for row in ddi_ids]

    random.shuffle(ddi_ids)
    total_ddis = len(ddi_ids)# Split and write DDI IDs to files
    for output_file, fraction in output_file_fraction_dict.items():
        count = int(total_ddis * fraction)
        ddi_subset = ddi_ids[:count]
        ddi_ids = ddi_ids[count:]

        with open(output_file, 'w') as f:
            f.write("ddi_id\\n")
            for ddi_id in ddi_subset:
                f.write(f"{ddi_id}\\n")

    conn.close()

    import sys as _sys
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {_sys.version.split()[0]}\\n")
    """
}