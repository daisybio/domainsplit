process RANDOM_DENOISE_SPLIT {
    tag "denoise ${change_fraction}${invert ? ' inverted' : ''}"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path 'cobinet.sqlite3'
    val change_fraction
    val invert

    output:
    path('*.csv'), emit: split_ddi_id_files
    val output_file_splits, emit: split_fractions
    path "versions.yml", emit: versions

    script:
    output_file_splits = []

    ["train", "optimization", "test"].each { name ->
        output_file_splits << ["${name}.csv", name]
    }

    """
    #!/usr/bin/env python3

    import sqlite3
    import random
    import os

    input_db_path = "cobinet.sqlite3"

    # Connect to the database
    conn = sqlite3.connect(input_db_path)
    # Fetch all domain-domain interaction IDs
    ddi_ids = conn.execute("SELECT id FROM domain_domain_interaction;")
    ddi_ids = [row[0] for row in ddi_ids]
    random.shuffle(ddi_ids)

    ddi_ids_positive = conn.execute("SELECT id FROM domain_domain_interaction WHERE negative = 0;")
    ddi_ids_positive = [row[0] for row in ddi_ids_positive]
    random.shuffle(ddi_ids_positive)

    ddi_ids_negative = conn.execute("SELECT id FROM domain_domain_interaction WHERE negative = 1;")
    ddi_ids_negative = [row[0] for row in ddi_ids_negative]
    random.shuffle(ddi_ids_negative)

    if ${invert ? 1 : 0}:
        print("Inverting denoising logic: changing positives to negatives.")
        replacement_value = 1
        ddi_ids_negative, ddi_ids_positive = ddi_ids_positive, ddi_ids_negative
    else:
        print("Standard denoising logic: changing negatives to positives.")
        replacement_value = 0


    change_count = int(len(ddi_ids_negative) * ${change_fraction})
    optimization_changed_ddis = set(ddi_ids_negative[:change_count])
    test_changed_ddis = set(ddi_ids[-change_count:])

    optimization_relevant_set = optimization_changed_ddis | set(ddi_ids_positive[:change_count])
    test_relevant_set = test_changed_ddis | set(ddi_ids_positive[-change_count:])

    train_set_unchanged = set(ddi_ids) - optimization_changed_ddis - test_changed_ddis

    print("Train set unchanged size:", len(train_set_unchanged))
    print("Optimization changed DDIs size:", len(optimization_changed_ddis))
    print("Test changed DDIs size:", len(test_changed_ddis))
    print("Optimization relevant set size:", len(optimization_relevant_set))
    print("Test relevant set size:", len(test_relevant_set))

    with open("train.csv", 'w') as f:
        f.write("ddi_id,override_negative\\n")
        for ddi_id in ddi_ids:
            if ddi_id in test_changed_ddis or ddi_id in optimization_changed_ddis:
                f.write(f"{ddi_id},{replacement_value}\\n")
            else:
                f.write(f"{ddi_id},\\n")

    with open("optimization.csv", 'w') as f:
        f.write("ddi_id,override_negative,is_evaluation_relevant\\n")
        for ddi_id in ddi_ids:
            if ddi_id in test_changed_ddis:
                f.write(f"{ddi_id},{replacement_value},0\\n")
            elif ddi_id in optimization_relevant_set:
                f.write(f"{ddi_id},,1\\n")
            else:
                f.write(f"{ddi_id},,\\n")

    with open("test.csv", 'w') as f:
        f.write("ddi_id,override_negative,is_evaluation_relevant\\n")
        for ddi_id in ddi_ids:
            if ddi_id in test_relevant_set:
                f.write(f"{ddi_id},,1\\n")
            else:
                f.write(f"{ddi_id},,\\n")


    conn.close()

    import sys as _sys
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {_sys.version.split()[0]}\\n")
    """

}