process DOWNLOAD_PFAM_ALIGNMENTS_BATCH {
    tag { "batch_${pfam_ids_list.size()}" }
    label 'process_low'
    conda "${moduleDir}/environment.yml"
    container "docker://konstantinpelz/domainsplit-general:1.0.0"

    maxRetries 3
    errorStrategy { task.attempt <= 3 ? 'retry' : 'ignore' }

    input:
    val pfam_ids_list

    output:
    path "*.alignment.full.gz", emit: alignments
    path "versions.yml", emit: versions

    script:
    def ids_joined = pfam_ids_list.collect { it.toString().strip() }.join('\n')
    def url_template = params.url_pfam_template
    """
#!/usr/bin/env python3
import sys, time, urllib.request, urllib.error, ssl
from concurrent.futures import ThreadPoolExecutor, as_completed

URL_TEMPLATE = "${url_template}"
PFAM_IDS = '''${ids_joined}'''.strip().split("\\n")

ctx = ssl.create_default_context()

def download_one(pfam_id):
    pfam_id = pfam_id.strip()
    if not pfam_id:
        return
    url = URL_TEMPLATE.replace("{pfam_id}", pfam_id)
    out_file = f"{pfam_id}.alignment.full.gz"
    req = urllib.request.Request(url, headers={"User-Agent": "domainsplit/1.0"})
    last_err = None
    for attempt in range(1, 6):
        try:
            with urllib.request.urlopen(req, timeout=120, context=ctx) as r, \\
                 open(out_file, "wb") as out:
                while True:
                    chunk = r.read(65536)
                    if not chunk:
                        break
                    out.write(chunk)
            return pfam_id
        except Exception as e:
            last_err = e
            time.sleep(min(2 ** attempt, 30))
    print(f"all retries failed for {pfam_id}: {last_err}", file=sys.stderr)
    return None

failed = []
with ThreadPoolExecutor(max_workers=8) as pool:
    futures = {pool.submit(download_one, pid): pid for pid in PFAM_IDS}
    for i, future in enumerate(as_completed(futures), 1):
        try:
            result = future.result()
        except Exception as e:
            result = None
            print(f"worker exception for {futures[future]}: {e}", file=sys.stderr)
        if result is None:
            failed.append(futures[future])
        if i % 20 == 0:
            print(f"  downloaded {i}/{len(PFAM_IDS)}", file=sys.stderr)

if failed:
    print(f"WARNING: {len(failed)}/{len(PFAM_IDS)} downloads failed: {failed[:10]}", file=sys.stderr)
if len(failed) == len(PFAM_IDS):
    print("ERROR: all downloads failed", file=sys.stderr)
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
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

def load_uniprot_id_mapping(path):
    mapping = {}
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            uniprot_id, id_type, symbol_name = line.split("\\t")
            if id_type.strip() == "UniProtKB-ID":
                mapping[symbol_name.strip()] = uniprot_id.strip()
    return mapping

def parse_stockholm_lightweight(fh):
    records = []
    for line in fh:
        if line.startswith('#') or line.startswith('//') or not line.strip():
            continue
        parts = line.split()
        if len(parts) != 2:
            continue
        name_field, raw_seq = parts
        if '/' not in name_field:
            continue
        name, coords = name_field.split('/', 1)
        if '-' not in coords:
            continue
        start, end = coords.split('-', 1)
        sequence = raw_seq.replace('-', '').replace('.', '').upper()
        if sequence:
            records.append((name, start, end, sequence))
    return records

def process_alignment(aln_path_str, id_map):
    aln = Path(aln_path_str)
    rows = []
    try:
        domain = aln.stem.split('.')[0]
        with gzip.open(aln, 'rt') as aln_fh:
            for name, start, end, sequence in parse_stockholm_lightweight(aln_fh):
                uniprot_id = id_map.get(name)
                if uniprot_id:
                    rows.append(f"{domain},{uniprot_id},{start},{end},{sequence}")
    except Exception:
        print(f"Warning: could not process alignment file {aln}", file=sys.stderr)
    return rows

print("Loading UniProt ID mapping...", file=sys.stderr)
id_map = load_uniprot_id_mapping('${uniprot_map_file.name}')

aln_files = sorted(Path('alignment_files').glob('*'))
print(f"Processing {len(aln_files)} alignment files with parallel workers...", file=sys.stderr)

num_workers = max(1, int(os.environ.get('NSLOTS', ${task.cpus})) - 1)
all_rows = []
with ProcessPoolExecutor(max_workers=num_workers) as pool:
    futures = {pool.submit(process_alignment, str(f), id_map): f for f in aln_files}
    done = 0
    for future in as_completed(futures):
        all_rows.extend(future.result())
        done += 1
        if done % 500 == 0:
            print(f"  processed {done}/{len(aln_files)} files", file=sys.stderr)

print(f"Writing {len(all_rows)} rows...", file=sys.stderr)
with gzip.open("$out_path", 'wt') as out:
    out.write("pfam_id,uniprot_id,start_pos,end_pos,sequence\\n")
    for row in all_rows:
        out.write(row + "\\n")

with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f"    python: {sys.version.split()[0]}\\n")
    """
}
