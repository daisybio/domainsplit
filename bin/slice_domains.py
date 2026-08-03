#! /usr/bin/env python3

import os
import shutil
import sys
from pathlib import Path
import argparse as ap

import Bio

import utils_struct


SOURCES = ("AF3", "RF2")


def parse_args():
    p = ap.ArgumentParser(description=__doc__, formatter_class=ap.RawDescriptionHelpFormatter)
    p.add_argument("--db_in", required=True)
    p.add_argument("--db_out", required=True)
    p.add_argument("--pdb_dir_af", required=True)
    p.add_argument("--pdb_dir_rf", required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process_name", required=True)
    return p.parse_args()



def slice_complex(conn, pdb_file: Path, source: str):
    """
    Slice every domain (single) and every matching DDI pair out of one
    predicted complex.
    """
    protein_id_a, protein_id_b = pdb_file.stem.split("_")
    chain_id_a, chain_id_b = "A", "B"

    domains_a = utils_struct.get_domain_mapping(conn, int(protein_id_a))
    domains_b = utils_struct.get_domain_mapping(conn, int(protein_id_b))

    if not domains_a or not domains_b:
        print(f"  WARNING: missing domain mapping for "
              f"{protein_id_a}/{protein_id_b} ({source}) -- "
              f"domains_a={len(domains_a)}, domains_b={len(domains_b)}", flush=True)


    
    for (domain_id_a, start_a, end_a) in domains_a:
        for (domain_id_b, start_b, end_b) in domains_b:
            ddi_id = utils_struct.check_ddi_exists(conn, domain_id_a, domain_id_b)
            if not ddi_id:
                continue  # this domain combination isn't a known DDI

            try:
                pdb_gz = utils_struct.ddi_pair_to_bytes(
                    str(pdb_file),
                    chain_id_a, start_a, end_a,
                    chain_id_b, start_b, end_b,
                )
            except Exception as exc:
                print(f"  WARNING: could not slice DDI pair "
                      f"{domain_id_a}/{domain_id_b} from {pdb_file}: {exc}", flush=True)
                continue

            utils_struct.store_domain_slice(conn, ddi_id, int(domain_id_a), int(domain_id_b), int(protein_id_a), int(protein_id_b), pdb_gz, source)


def _write_versions(versions_path: str, process_name: str) -> None:
    with open(versions_path, "w") as fh:
        fh.write(f'"{process_name}":\n')
        fh.write(f"    python: {sys.version.split()[0]}\n")
        fh.write(f"    biopython: {Bio.__version__}\n")


def enrich_db_with_predictions(conn, pdb_dir: Path, source: str):
    
    pdb_files = list(pdb_dir.glob(f"*.pdb"))
    print(f"[slice_domains] source={source}: {len(pdb_files)} predicted "
          f"complexes to slice", flush=True)

    for pdb_path in pdb_files:
        slice_complex(conn, pdb_path, source)     
        conn.commit()
   


def main():
    args = parse_args()
    
    pdb_dir_af = Path(args.pdb_dir_af) # contains predicted PDBs for AF3
    pdb_dir_rf = Path(args.pdb_dir_rf) # contains predicted PDBs for RF

    shutil.copy(args.db_in, args.db_out)
    conn = utils_struct.connect_db(args.db_out)

    # Check if output directories exist, if not skip slicing for that source
    if not pdb_dir_af.exists():
        print(f"[slice_domains] WARNING: AF3 PDB directory {pdb_dir_af} does not exist, skipping AF3 slicing", flush=True)
    else:
        enrich_db_with_predictions(conn, pdb_dir_af, "AF3")
    if not pdb_dir_rf.exists():
        print(f"[slice_domains] WARNING: RF PDB directory {pdb_dir_rf} does not exist, skipping RF slicing", flush=True)
    else:
        enrich_db_with_predictions(conn, pdb_dir_rf, "RF")
    conn.close()

    _write_versions(args.versions, args.process_name)





if __name__ == "__main__":
    main()