process PROTT5_EMBEDDINGS_CHUNK {
    tag { tag_str }
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    val(protein_ids)

    output:
    path "embeddings_chunk_${tag_str}.h5.gz", emit: embeddings_chunk
    path "versions.yml", emit: versions

    script:
    tag_str = "${protein_ids[0]}-${protein_ids[-1]}"
    query_string = protein_ids.collect { "(accession:$it)" }.join("OR")

    """
    #!/usr/bin/env python3

    import requests
    import time

    output_path = "embeddings_chunk_${tag_str}.h5.gz"
    query_string = "${query_string}"
    params = {
        "format": "h5",
        "query": query_string,
    }

    print("Submitting request to UniProt for chunk: ${tag_str}")
    url = "https://rest.uniprot.org/uniprotkb/download/run"
    response = requests.post(url, params=params)

    try:
        response.raise_for_status()
    except requests.HTTPError as e:
        print(response.status_code)
        print(response.text)
        print("Failed to submit request to UniProt:", e)
    else:
        job_id = response.json()["jobId"]

        def get_job_status():
            status_url = f"https://rest.uniprot.org/uniprotkb/download/status/{job_id}"
            r = requests.get(status_url)
            r.raise_for_status()
            status = r.json()["jobStatus"]
            print(f"Job {job_id} status: {status}")
            return status

        # wait until the request is done
        time.sleep(2)
        while get_job_status() == "RUNNING":
            time.sleep(0.2)

        # download the file
        print("Downloading results for job:", job_id)
        results_url = f"https://rest.uniprot.org/uniprotkb/download/results/{job_id}"
        with requests.get(results_url, stream=True) as r, open(output_path, "wb") as out_file:
            r.raise_for_status()
            for chunk in r.iter_content(chunk_size=65536):
                out_file.write(chunk)

    import sys as _sys
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {_sys.version.split()[0]}\\n")
        f.write(f"    requests: {requests.__version__}\\n")
    """
}

process JOIN_HDF_FILES {
    tag "join_prott5_chunks"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    path "chunk*"

    output:
    path "embeddings.h5", emit: embeddings
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    from pathlib import Path
    import h5py
    import gzip

    with h5py.File("embeddings.h5", "w") as out_h5:
        for chunk_file in Path().glob("chunk*"):
            print("Combining chunk file:", chunk_file)
            # open chunk file with gzip decompression
            with gzip.open(chunk_file, "rb") as gz_file:
                with h5py.File(gz_file, "r") as in_h5:
                    for key in in_h5.keys():
                        in_h5.copy(key, out_h5)

    import sys as _sys
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {_sys.version.split()[0]}\\n")
        f.write(f"    h5py: {h5py.__version__}\\n")
    """
}

process DOWNLOAD_PROTT5_EMBEDDINGS_COMPLETE {
    tag "prott5_complete"
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    maxRetries 3
    errorStrategy "ignore"

    output:
    path(output_path), emit: embeddings
    path "versions.yml", emit: versions

    script:
    download_url = params.url_prott5_embeddings_complete
    output_path = "uniprot_sprot_per_residue.h5.gz"
    """
    wget -O "${output_path}" "${download_url}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | head -n1 | awk '{print \$3}')
    END_VERSIONS
    """
}

// Helper to retrieve ProtT5 per-residue embeddings.
//
// Two paths exist:
//   1. `download_prott5_complete = true` (default): grab the full UniProt-SwissProt per-residue HDF5
//      from EBI (`params.url_prott5_embeddings_complete`). Single big HTTP GET; reliable; large file.
//   2. `download_prott5_complete = false`: post chunked queries to UniProt's async REST API for only
//      the required protein IDs, then join chunked HDF5s. Bandwidth-friendly when the relevant
//      protein subset is small, but the chunked API has historically been slow/unstable.
//
// Default stays on the complete download to preserve known-good behaviour. Flip when you want
// to filter at fetch time (e.g. species-specific runs).
def download_prott5_embeddings(protein_ids) {
    if (params.download_prott5_complete) {
        return DOWNLOAD_PROTT5_EMBEDDINGS_COMPLETE().embeddings
    }
    def protein_id_sublists = protein_ids.buffer(size: 100, remainder: true)
    def chunks = PROTT5_EMBEDDINGS_CHUNK(protein_id_sublists)
    def embeddings = JOIN_HDF_FILES(chunks.embeddings_chunk.collect())
    return embeddings.embeddings
}
