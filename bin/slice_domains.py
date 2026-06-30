#! /usr/bin/env python3

import shutil
import sys
from pathlib import Path
import argparse as ap

import Bio

import utils_struct


SOURCES = ("AF3", "RF2")


def parse_args():
    p = ap.ArgumentParser(description=__doc__, formatter_class=ap.RawDescriptionHelpFormatter)
    p.add_argument("--db_in", required=True, help="SQLite DB with complex_chain_map populated by the predictors")
    p.add_argument("--db_out", required=True)
    p.add_argument("--pdb_dir", required=True)
    p.add_argument("--versions", required=True)
    p.add_argument("--process_name", required=True)
    return p.parse_args()



def slice_complex(conn, pdb_dir: Path, source: str, protein_id_a: int, uniprot_id_a: str, chain_id_a: str, protein_id_b: int, uniprot_id_b: str, chain_id_b: str):
    """
    Slice every domain (single) and every matching DDI pair out of one
    predicted complex.  Returns (n_domains_stored, n_ddis_stored).
    """
    pdb_path = pdb_dir / f"{uniprot_id_a}_{uniprot_id_b}_{source}.pdb"
    if not pdb_path.exists():
        print(f"  WARNING: missing predicted PDB for {uniprot_id_a}/{uniprot_id_b} "
              f"({source}): {pdb_path}", flush=True)
        return 0, 0

    domains_a = utils_struct.get_domain_mapping(conn, protein_id_a)
    domains_b = utils_struct.get_domain_mapping(conn, protein_id_b)

    if not domains_a or not domains_b:
        print(f"  WARNING: missing domain mapping for "
              f"{uniprot_id_a}/{uniprot_id_b} ({source}) -- "
              f"domains_a={len(domains_a)}, domains_b={len(domains_b)}", flush=True)


    n_ddis = 0
    for (domain_id_a, start_a, end_a) in domains_a:
        for (domain_id_b, start_b, end_b) in domains_b:
            ddi_id = utils_struct.check_ddi_exists(conn, domain_id_a, domain_id_b)
            if not ddi_id:
                continue  # this domain combination isn't a known DDI

            try:
                pdb_gz = utils_struct.ddi_pair_to_bytes(
                    str(pdb_path),
                    chain_id_a, start_a, end_a,
                    chain_id_b, start_b, end_b,
                )
            except Exception as exc:
                print(f"  WARNING: could not slice DDI pair "
                      f"{domain_id_a}/{domain_id_b} from {pdb_path}: {exc}", flush=True)
                continue

            utils_struct.store_domain_slice(conn, protein_id_a, protein_id_b, ddi_id, pdb_gz, source)
            n_ddis += 1


def _write_versions(versions_path: str, process_name: str) -> None:
    with open(versions_path, "w") as fh:
        fh.write(f'"{process_name}":\n')
        fh.write(f"    python: {sys.version.split()[0]}\n")
        fh.write(f"    biopython: {Bio.__version__}\n")


def main():
    args = parse_args()
    pdb_dir = Path(args.pdb_dir) # contains predicted PDBs

    shutil.copy(args.db_in, args.db_out)
    conn = utils_struct.connect_db(args.db_out)
    utils_struct.create_tables(conn)

    total_domains = total_ddis = total_complexes = 0

    for source in SOURCES:
        complexes = utils_struct.get_predicted_complexes(conn, source)
        print(f"[slice_domains] source={source}: {len(complexes)} predicted "
              f"complexes to slice", flush=True)

        for (protein_id_a, uniprot_id_a, chain_a, protein_id_b, uniprot_id_b, chain_b) in complexes:
            slice_complex(
                conn, pdb_dir, source,
                protein_id_a, uniprot_id_a, chain_a,
                protein_id_b, uniprot_id_b, chain_b,
            )
            conn.commit()

    conn.close()
    _write_versions(args.versions, args.process_name)
    print(f"[slice_domains] done: {total_complexes} complexes processed, "
          f"{total_domains} domain slices, {total_ddis} DDI-pair slices stored",
          flush=True)





if __name__ == "__main__":
    main()