process DOWNLOAD_PFAM_ALIGNMENT {
    tag { pfam_id }
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    maxRetries 3
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

    input:
    val pfam_id

    output:
    path "${pfam_id}.alignment.full.gz", emit: alignment
    path "versions.yml", emit: versions

    script:
    pfam_id = pfam_id.strip()
    def download_url = params.url_pfam_template.replace("{pfam_id}", pfam_id)
    """
    #!/usr/bin/env python3
    import sys, time, urllib.request, urllib.error, ssl

    OUTPUT_FILE = "${pfam_id}.alignment.full.gz"
    DOWNLOAD_URL = "${download_url}"

    ctx = ssl.create_default_context()
    req = urllib.request.Request(DOWNLOAD_URL, headers={"User-Agent": "domainsplit/1.0"})

    last_err = None
    for attempt in range(1, 6):
        try:
            with urllib.request.urlopen(req, timeout=120, context=ctx) as r, open(OUTPUT_FILE, "wb") as out:
                while True:
                    chunk = r.read(65536)
                    if not chunk:
                        break
                    out.write(chunk)
            break
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            last_err = e
            print(f"attempt {attempt} failed for ${pfam_id}: {e}", file=sys.stderr)
            time.sleep(2 ** attempt)
    else:
        print(f"all retries failed for ${pfam_id}: {last_err}", file=sys.stderr)
        sys.exit(1)

    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
    """
}

process CREATE_PROTEIN_DOMAIN_MAPPING {
    tag { "${uniprot_map_file.simpleName}" }
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    input:
    path uniprot_map_file
    path "alignment_files/*"

    output:
    path out_path, emit: mapping
    path "versions.yml", emit: versions

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

    import Bio as _Bio
    with open("versions.yml", "w") as f:
        f.write('"${task.process}":\\n')
        f.write(f"    python: {sys.version.split()[0]}\\n")
        f.write(f"    biopython: {_Bio.__version__}\\n")
    """
}