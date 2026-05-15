process PROTT5_EMBEDDINGS_CHUNK {
    tag { tag_str }
    conda "${moduleDir}/environment.yml"
    maxForks 20

    input:
    val(protein_ids)

    output:
    path "embeddings_chunk_${tag_str}.h5.gz", emit: embeddings_chunk

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
    """
}

process JOIN_HDF_FILES {
    conda "${moduleDir}/environment.yml"

    input:
    path "chunk*"

    output:
    path "embeddings.h5"

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
    """
}

process DOWNLOAD_PROTT5_EMBEDDINGS_COMPLETE {
    conda "${moduleDir}/environment.yml"

    maxRetries 3
    errorStrategy "ignore"

    output:
    path(output_path)

    script:
    download_url = "https://ftp.ebi.ac.uk/pub/contrib/UniProt/embeddings/current_release/uniprot_sprot/per-residue.h5"
    output_path = "uniprot_sprot_per_residue.h5.gz"
    """
    wget -O "${output_path}" "${download_url}"
    """
}

def download_prott5_embeddings(protein_ids) {
    if(true)
        //return file("https://ftp.ebi.ac.uk/pub/contrib/UniProt/embeddings/current_release/uniprot_sprot/per-residue.h5")
        return DOWNLOAD_PROTT5_EMBEDDINGS_COMPLETE()
    protein_id_sublists = protein_ids.buffer(size:100, remainder:true)
    chunks = PROTT5_EMBEDDINGS_CHUNK(protein_id_sublists)
    embeddings = JOIN_HDF_FILES(chunks.collect())
    return embeddings
}
