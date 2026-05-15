process DOWNLOAD_PFAM_ALIGNMENT {
    conda "${moduleDir}/environment.yml"
    tag { pfam_id }

    maxRetries 3
    maxForks 100
    errorStrategy "ignore"

    input:
    val pfam_id

    output:
    path "${pfam_id}.alignment.full.gz"

    script:
    pfam_id = pfam_id.strip()
    def download_url = params.url_pfam_template.replace("{pfam_id}", pfam_id)
    """
    OUTPUT_FILE='${pfam_id}.alignment.full.gz'
    DOWNLOAD_URL='$download_url'

    wget -O "\$OUTPUT_FILE" "\$DOWNLOAD_URL"
    """
}

process CREATE_PROTEIN_DOMAIN_MAPPING {
    conda "${moduleDir}/environment.yml"

    input:
    path uniprot_map_file
    path "alignment_files/*"

    output:
    path out_path

    script:
    out_path = 'protein_domain_mapping.csv.gz'
    """
    #!/usr/bin/env python3
    import gzip, sys, os
    from Bio import AlignIO
    from pathlib import Path

    def load_uniprot_id_mapping(path):
        mapping = {}
        with gzip.open(path, 'rt') as fh:
            for line in fh:
                uniprot_id, id_type, symbol_name = line.split("\t")
                if id_type.strip() == "UniProtKB-ID":
                    mapping[symbol_name.strip()] = uniprot_id.strip()
        return mapping

    print("Loading UniProt ID mapping...", file=sys.stderr)
    id_map = load_uniprot_id_mapping('${uniprot_map_file.name}')

    print("Processing alignment files...", file=sys.stderr)
    with gzip.open("$out_path", 'wt') as out:
        out.write("pfam_id,uniprot_id,start_pos,end_pos,sequence\\n")
        for aln in Path('alignment_files').glob('*'):
            try:
                with gzip.open(aln, 'rt') as aln_fh:
                    alignment = AlignIO.read(aln_fh, 'stockholm')
                    domain = aln.stem.split('.')[0]
                    for record in alignment:
                        uniprot_name = record.id.split('/')[0]
                        uniprot_id = id_map.get(uniprot_name)
                        start, end = record.id.split("/")[1].split("-")
                        sequence = str(record.seq.replace("-","")).upper()
                        if uniprot_id and (domain, uniprot_id):
                            out.write(f"{domain},{uniprot_id},{start},{end},{sequence}\\n")
            except Exception:
                print(f"Warning: could not process alignment file {aln}", file=sys.stderr)
    """
}