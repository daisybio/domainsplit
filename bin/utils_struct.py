#!/usr/bin/env python3

import gzip
import io
import sqlite3

from Bio.PDB.PDBIO import Select
from Bio.PDB.SASA import ShrakeRupley


# ---------------------------------------------------------------------------
# Instance-id / structures.h5 key helpers (shared by enrich_structural_af.py,
# build_scoring_matrix.py, score_ddi.py)
# ---------------------------------------------------------------------------

def extract_data_from_instance_id(instance_id: str):
    """instance_id is pfam_id + '_' + uniprot_id + '_' + start + '_' + end."""
    parts = instance_id.split("_")
    if len(parts) != 4:
        raise ValueError(f"instance_id {instance_id!r} does not have 4 parts separated by '_'")
    pfam_id, uniprot_id, start, end = parts
    return pfam_id, uniprot_id, int(start), int(end)


def structure_h5_keys(instance_id_a: str, instance_id_b: str):
    """(group_key, dataset_key) into structures.h5 for one instance pair."""
    pfam_a, *_ = extract_data_from_instance_id(instance_id_a)
    pfam_b, *_ = extract_data_from_instance_id(instance_id_b)
    return f"{pfam_a}_{pfam_b}", f"{instance_id_a}__{instance_id_b}"


def get_structure_bytes(h5file, instance_id_a: str, instance_id_b: str):
    """Gzip-compressed PDB bytes for one instance pair, or None if absent."""
    group_key, dset_key = structure_h5_keys(instance_id_a, instance_id_b)
    grp = h5file.get(group_key)
    if grp is None or dset_key not in grp:
        return None
    return grp[dset_key][()].tobytes()



# ---------------------------------------------------------------------------
# Structure Helpers
# ---------------------------------------------------------------------------


def bytes_to_stringio(blob: bytes) -> io.StringIO:
    """
    Decompress gzip-compressed PDB bytes into an in-memory text stream.
    PDBParser.get_structure() accepts a file-like object directly.
    """
    return io.StringIO(gzip.decompress(blob).decode("utf-8"))


class DdiPairSelect(Select):

    def __init__(self, chain_id_a: str, start_a: int, end_a: int, chain_id_b: str, start_b: int, end_b: int) -> None:
        self.window_a = (chain_id_a, start_a, end_a)
        self.window_b = (chain_id_b, start_b, end_b)
 
    def _in_window(self, residue, window) -> bool:
        chain_id, start, end = window
        return (
            residue.get_parent().id == chain_id
            and start <= residue.id[1] <= end
        )
 
    def accept_residue(self, residue):
        if self._in_window(residue, self.window_a):
            return 1
        if self._in_window(residue, self.window_b):
            return 1
        return 0
 


#---------------------------------------------------------------------------
# Scoring Helpers
#---------------------------------------------------------------------------


AA_3 = [
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
]
AA_SET = set(AA_3)

BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}

NO_HB  = 3.5   # H-bond:     N-O distance <= 3.5 Å
NO_SB  = 5.5   # Salt bridge: N-O distance <= 5.5 Å
CC_VDW = 5.0   # vdW:        C-C distance <= 5.0 Å



def calculate_sasa_residue_level(domain, chain_id):
    sr = ShrakeRupley()
    if chain_id:
        domain = domain[0][chain_id]  # Get the specific chain from the structure
    sr.compute(domain, level="R")  # Compute SASA at the residue level
    sasa_values = {}
    for residue in domain.get_residues():
        sasa_values[residue.get_id()] = residue.sasa
    return sasa_values


def calculate_rsa_residue_level(domain, chain_id):
    sasa_residue = calculate_sasa_residue_level(domain, chain_id)

    # MAxSASA values by Tien et al. 2013, "Maximum allowed solvent accessibilities of residues in proteins" (https://doi.org/10.1002/prot.24286)
    max_sasa_values = {
        'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0, 'CYS': 167.0,
        'GLN': 223.0, 'GLU': 225.0, 'GLY': 104.0, 'HIS': 224.0, 'ILE': 197.0,
        'LEU': 201.0, 'LYS': 236.0, 'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0,
        'SER': 155.0, 'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0
    }

    rsa_residue = {}
    for residue in domain.get_residues():
        resname = residue.get_resname()
        rid = residue.get_id()
        sasa_value = sasa_residue.get(rid, 0)
        max_sasa = max_sasa_values.get(resname, None)
        if max_sasa is not None and max_sasa > 0:
            rsa_residue[rid] = sasa_value / max_sasa
        else:
            rsa_residue[rid] = None

    return rsa_residue



def extract_domain_residues(structure, chain_id):
    """
    Return standard amino-acid residue objects from chain_id within
    the structure. Non-standard residues are ignored.
    """
    try:
        chain = structure[0][chain_id]
    except (KeyError, StopIteration, IndexError):
        return []
    return [
        r for r in chain.get_residues()
        if r.get_resname() in AA_SET
    ]



# ---------------------------------------------------------------------------
# domain_structure access -- rewritten for the H5-backed schema (no pdb_gz,
# no source column: RF is gone, bytes live in structures.h5). Both exclude
# is_mock rows: see enrich_structural_af.py's mock-fallback path.
# ---------------------------------------------------------------------------

def get_positive_train_structures(db_path, method):
    """Real (non-mock) structures for POSITIVE DDIs in this method's TRAIN
    split only -- feeds this method's own potential in build_scoring_matrix.py.
    Must never include an instance outside this method's train partition, or
    a test-time DDI scored later under this same matrix would be evaluated
    against a background that had already seen it.
    Returns [(ddi_id, instance_id_a, instance_id_b)]."""
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA synchronous  = OFF")
    rows = conn.execute("""
        SELECT DISTINCT m.ddi_id, m.instance_id_a, m.instance_id_b
        FROM ddi_split_membership m
        JOIN domain_domain_interaction ddi ON ddi.id = m.ddi_id
        WHERE m.method = ? AND m.split = 'train' AND ddi.negative = 0 AND m.is_mock = 0
    """, (method,)).fetchall()
    conn.close()
    return rows


def get_structures_for_method(db_path, method):
    """Real (non-mock) structures for every DDI this method places in ANY
    split -- feeds score_ddi.py, which scores the whole method.
    Returns [(ddi_id, instance_id_a, instance_id_b, split)]."""
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA synchronous  = OFF")
    rows = conn.execute("""
        SELECT DISTINCT m.ddi_id, m.instance_id_a, m.instance_id_b, m.split
        FROM ddi_split_membership m
        WHERE m.method = ? AND m.is_mock = 0
    """, (method,)).fetchall()
    conn.close()
    return rows

