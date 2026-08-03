process PREDICT_COMPLEX_AF {
    tag "predict_complex_af"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path dbstruct, stageAs: 'input.dbstruct.sqlite3'
    path meta_ppis

    output:
    path "pdb_files_af/", emit: pdb_files
    path "versions.yml",  emit: versions

    script:
    """
    #!/usr/bin/env python3
    import glob
    import os
    import sys
    import sqlite3
    import pandas as pd
    import zstandard
    import Bio
    from Bio.PDB import MMCIFParser, PDBIO

    OUTDIR = "pdb_files_af"
    os.makedirs(OUTDIR, exist_ok=True)

    def convert_cif_to_pdb(cif_path: str, pdb_path: str) -> None:
        try:
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("model", cif_path)
            io_obj = PDBIO()
            io_obj.set_structure(structure)
            io_obj.save(pdb_path)
        except Exception as e:
            print(f"Error converting {cif_path} to {pdb_path}: {e}", flush=True)

    def decompress_zst(src_path: str, dst_path: str) -> None:
        dctx = zstandard.ZstdDecompressor()
        with open(src_path, "rb") as fh_in, open(dst_path, "wb") as fh_out:
            dctx.copy_stream(fh_in, fh_out)

    # Load metadata
    meta = pd.read_csv("${meta_ppis}")
    meta_index = {}
    for _, row in meta.iterrows():
        key_fwd = (row["uniprot_id_a"], row["uniprot_id_b"])
        key_rev = (row["uniprot_id_b"], row["uniprot_id_a"])
        path_to_file = row["path"]
        has_af_model = bool(row["has_af_model"])
        order = bool(row["ordering"])
        # need to preserve ordering information to know which protein is A and which is B
        # If the order is true and we have forward key, we set true
        # if the order is true and we have reverse key, we set false
        # if the order is false and we have forward key, we set false
        # if the order is false and we have reverse key, we set true
        if order:
            meta_index[key_fwd] = (path_to_file, has_af_model, True)
            meta_index[key_rev] = (path_to_file, has_af_model, False)
        else:
            meta_index[key_fwd] = (path_to_file, has_af_model, False)
            meta_index[key_rev] = (path_to_file, has_af_model, True)

        # meta_index[key_fwd] = (path_to_file, has_af_model, )
        # meta_index[key_rev] = (path_to_file, has_af_model, )

    # Get PPIs in database, join to get uniprot ids
    con = sqlite3.connect(f"file:input.dbstruct.sqlite3?mode=ro", uri=True)
    ppi_rows = con.execute(
        "SELECT p1.uniprot_id, p2.uniprot_id, ppi.protein_id_a, ppi.protein_id_b FROM protein_protein_interaction ppi "
        "JOIN protein p1 ON ppi.protein_id_a = p1.id "
        "JOIN protein p2 ON ppi.protein_id_b = p2.id "
        "WHERE p1.uniprot_id IS NOT NULL AND p2.uniprot_id IS NOT NULL"
    ).fetchall()

    con.close()

    print(f"predict_complex_af: {len(ppi_rows)} PPIs to resolve", flush=True)

    n_ok = 0
    n_missing_meta = 0
    n_missing_file = 0
    n_ambiguous = 0

    for uniprot_id_a, uniprot_id_b, protein_id_a, protein_id_b in ppi_rows:
        path, has_af_model, order = meta_index.get((uniprot_id_a, uniprot_id_b), (None, None, None))
        if path is None:
            print(f"WARNING: no metadata for pair ({uniprot_id_a}, {uniprot_id_b})", flush=True)
            n_missing_meta += 1
            continue

        # Path to output PDB file
        # Use protein ids for easier downstream mapping, avoids joining step
        if not order:
            # Ensure first protein is chain A
            protein_id_a, protein_id_b = protein_id_b, protein_id_a
        pdb_path = os.path.join(OUTDIR, f"{protein_id_a}_{protein_id_b}.pdb")

        if bool(has_af_model):
            candidates = glob.glob(os.path.join(path, "*.pdb.zst"))
        else:
            candidates = glob.glob(os.path.join(path, "*.cif"))

        if len(candidates) == 0:
            print(
                f"WARNING: no source file found in {path} "
                f"for ({uniprot_id_a}, {uniprot_id_b})", flush=True,
            )
            n_missing_file += 1
            continue
        if len(candidates) > 1:
            print(
                f"WARNING: {len(candidates)} candidate files in {path} "
                f"for ({uniprot_id_a}, {uniprot_id_b}), using first: {candidates[0]}",
                flush=True,
            )
            n_ambiguous += 1

        src_path = candidates[0]

        if bool(has_af_model):
            decompress_zst(src_path, pdb_path)
        else:
            convert_cif_to_pdb(src_path, pdb_path)

        n_ok += 1

    print(
        f"predict_complex_af: done -> ok={n_ok} missing_meta={n_missing_meta} "
        f"missing_file={n_missing_file} ambiguous={n_ambiguous}",
        flush=True,
    )

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\\n")
        f.write(f"    biopython: {Bio.__version__}\\n")
        f.write(f"    zstandard: {zstandard.__version__}\\n")
    """

    stub:
    """
    mkdir -p pdb_files_af
    echo '"${task.process}":' > versions.yml
    echo '    stub: "true"' >> versions.yml
    """
}