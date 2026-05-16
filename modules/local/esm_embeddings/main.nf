process FILTER_SEQUENCES {
    tag { "${protein_domain_map.simpleName}" }
    label 'process_medium'
    conda "${moduleDir}/environment.yml"

    input:
    path protein_domain_map
    path uniprotkb_database

    output:
    tuple val(protein_meta), path("uniprot_filtered.fasta.gz"), emit: protein_sequences
    tuple val(domain_meta), path("domain_sequences.fasta.gz"), emit: domain_sequences
    path "versions.yml", emit: versions

    script:
    protein_meta = [id: "protein_sequences"]
    domain_meta = [id: "domain_sequences"]
    """
    #!/usr/bin/env python3
    import gzip
    import pandas as pd
    from Bio import SeqIO

    with gzip.open("${protein_domain_map}","rt") as protein_domain_map_text:
        pd_map = pd.read_csv(protein_domain_map_text)

    with gzip.open("domain_sequences.fasta.gz", "wt") as domain_sequences_fasta:
        for row in pd_map.itertuples():
            domain_sequences_fasta.write(f">{row.pfam_id}_{row.uniprot_id}_{row.start_pos}_{row.end_pos}\\n{row.sequence}\\n")


    uniprot_ids = set(pd_map.uniprot_id)
    with gzip.open("uniprot_filtered.fasta.gz", "wt") as uniprot_filtered_fasta:
        with gzip.open("${uniprotkb_database}", "rt") as uniprotkb_fasta:
            for record in SeqIO.parse(uniprotkb_fasta, "fasta"):
                record.id = record.id.split("|")[1]
                if record.id in uniprot_ids:
                    record.description=""
                    SeqIO.write(record, uniprot_filtered_fasta, "fasta")

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        import sys, Bio, pandas
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    biopython: {Bio.__version__}\\n")
        f.write(f"    pandas: {pandas.__version__}\\n")
    """
}

process GENERATE_ESM_EMBEDDINGS {
    tag { meta.id }
    label 'process_gpu_large'
    maxForks 10
    conda "${moduleDir}/environment.yml"
    // NOTE: queue/clusterOptions/memory/time will move to conf/slurm.config in phase 5.
    // Kept inline for now so phase 1 is behavior-preserving.
    queue "shared-gpu"
    memory '20 GB'
    time 24.h
    clusterOptions "--gpus-per-task=1 --ntasks=1  --qos=limitgpus --nodelist=jlab-gpu01.exbio.wzw.tum.de,gpu02.exbio.wzw.tum.de,compms-gpu-2.exbio.wzw.tum.de"
    //  compms-gpu-1.exbio.wzw.tum.de is currently running having memory ECC issues, so we exclude it for now

    input:
    tuple val(meta), path(input_fasta)

    output:
    tuple val(meta), path("esm_embeddings.h5"), emit: esm_embeddings
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import os
    # HF_TOKEN must be provided by the executor environment (e.g. Nextflow secret or cluster env).
    # Source: see CLAUDE.md / params.hf_token_env.
    hf_token = os.environ.get("HF_TOKEN")
    if hf_token:
        os.system(f"hf auth login --token {hf_token} --add-to-git-credential")
    os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"

    from esm.models.esmc import ESMC
    from esm.models.esm3 import ESM3
    from esm.sdk.api import ESMProtein, LogitsConfig
    import h5py
    import sys
    from tqdm import tqdm
    from Bio import SeqIO
    import gzip
    import torch

    gpu = torch.device("cuda") if torch.cuda.is_available() else torch.device("cpu")
    #cpu = torch.device("cpu")
    config = LogitsConfig(sequence=True, return_embeddings=True)

    def create_embedding(client, sequence):
        with torch.inference_mode():
            protein = ESMProtein(sequence=sequence)
            encoded = client.encode(protein)
            embedding_tensor = client.logits(encoded, config).embeddings
            result = embedding_tensor.cpu().squeeze(dim=0) # remove batch dimension

            del embedding_tensor
            del encoded
            return result

    def yield_embeddings(client_gpu, client_cpu):
        # load and sort all records by sequence length (ascending)
        with gzip.open("$input_fasta", "rt") as input_fasta:
            records = list(SeqIO.parse(input_fasta, "fasta"))
        records.sort(key=lambda r: len(r.seq))

        for record in tqdm(records):
            seq = str(record.seq)
            for attempt in range(3): # try up to 3 times to create the embedding, in case of OOM errors
                try:
                    embedding = create_embedding(client_gpu, seq)
                    yield record.id, embedding
                    break # success, move to next sequence
                except torch.OutOfMemoryError:
                    print("CUDA OOM detected.")
                    print(f"Sequence length: {len(seq)}")
                    print(f"Attempt {attempt+1} of 3")
                    if attempt >= 2:
                        print(f"Failed to create embedding for sequence {record.id} after 3 attempts, skipping.")
                        return # give up after 3 attempts

    #def yield_embeddings(client):
    #    with gzip.open("$input_fasta", "rt") as input_fasta:
    #        for record in tqdm(SeqIO.parse(input_fasta, "fasta")):
    #            if len(record.seq) >= 5000:
    #                print(f"Skipping sequence {record.id}")
    #                continue
    #            yield record.id, create_embedding(client, str(record.seq))

    with h5py.File(f'esm_embeddings.h5', 'w') as esm_file:
        client_gpu = ESM3.from_pretrained("esm3-open", device=gpu)
        #client_cpu = ESM3.from_pretrained("esm3-open", device=cpu)
        client_cpu=None

        print("Generating ESM3 embeddings")
        for seq_name, embedding in yield_embeddings(client_gpu, client_cpu):
            esm_file.create_dataset(f"{seq_name}/esm3", data=embedding)

        client_gpu = ESMC.from_pretrained("esmc_600m", device=gpu)
        #client_cpu = ESMC.from_pretrained("esmc_600m", device=cpu)
        client_cpu = None

        print("Generating ESMC embeddings")
        for seq_name, embedding in yield_embeddings(client_gpu, client_cpu):
            esm_file.create_dataset(f"{seq_name}/esmc", data=embedding)

    import esm as _esm
    import torch as _torch
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    esm: {_esm.__version__}\\n")
        f.write(f"    torch: {_torch.__version__}\\n")
    """
}

process AVERAGE_POOL_EMBEDDINGS {
    tag { "${input_embeddings.simpleName}" }
    label 'process_low'
    conda "${moduleDir}/environment.yml"

    input:
    path("input_embeddings.h5")

    output:
    path("output_embeddings.h5"), emit: averaged_embeddings
    path "versions.yml", emit: versions

    script:
    """
    #!/usr/bin/env python3
    import h5py
    import numpy as np

    with h5py.File("input_embeddings.h5", "r") as input_h5, \
        h5py.File("output_embeddings.h5", "w") as output_h5:

        def copy_to_new_h5(path, item):
            if isinstance(item, h5py.Dataset):
                avg_embedding = np.mean(item, axis=0)
                output_h5.create_dataset(path, data=avg_embedding)

        input_h5.visititems(copy_to_new_h5)

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    h5py: {h5py.__version__}\\n")
        f.write(f"    numpy: {np.__version__}\\n")
    """
}


def generate_esm_embeddings(uniprotkb_database, protein_domain_map) {
    filter_result = FILTER_SEQUENCES(protein_domain_map, uniprotkb_database)
    protein_fasta = filter_result.protein_sequences
    domain_fasta = filter_result.domain_sequences

    esm_result = GENERATE_ESM_EMBEDDINGS(protein_fasta.mix(domain_fasta))

    protein_embeddings = esm_result.esm_embeddings.filter { it[0].id == "protein_sequences" }.map { it[1] }
    domain_embeddings = esm_result.esm_embeddings.filter { it[0].id == "domain_sequences" }.map { it[1] }

    return [protein_embeddings, AVERAGE_POOL_EMBEDDINGS(domain_embeddings).averaged_embeddings]
}
