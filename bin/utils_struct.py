#!/usr/bin/env python3

import gzip
import io
import os
import tempfile
import sqlite3

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO, Select
from Bio.PDB.SASA import ShrakeRupley


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


# ---------------------------------------------------------------------------
# Select domain coordinates class
# ---------------------------------------------------------------------------

class DomainSelect(Select):
    """Accept only residues in [start_residue, end_residue] on chain_id."""

    def __init__(self, chain_id: str, start_residue: int, end_residue: int):
        self.chain_id      = chain_id
        self.start_residue = start_residue
        self.end_residue   = end_residue

    def accept_residue(self, residue):
        if residue.get_parent().id == self.chain_id:
            if self.start_residue <= residue.id[1] <= self.end_residue:
                return 1
        return 0


# ---------------------------------------------------------------------------
# PDB Helpers
# ---------------------------------------------------------------------------

def extract_domain_pdb(pdb_file: str, chain_id: str, start_residue: int, end_residue: int, output_file: str) -> None:
    """Write domain slice to *output_file* in PDB-format"""
    parser    = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    io_obj    = PDBIO()
    io_obj.set_structure(structure)
    io_obj.save(output_file, DomainSelect(chain_id, start_residue, end_residue))



def domain_to_bytes(pdb_file: str, chain_id: str, start_residue: int, end_residue: int) -> bytes:
    """
    Extract the domain slice and return it as gzip-compressed PDB.
    -> stored in domain_structure.pdb_gz
    """
    parser    = PDBParser(QUIET=True)
    structure = parser.get_structure("domain", pdb_file)

    buf    = io.StringIO()
    io_obj = PDBIO()
    io_obj.set_structure(structure)
    io_obj.save(buf, DomainSelect(chain_id, start_residue, end_residue))

    raw = buf.getvalue().encode("utf-8")
    return gzip.compress(raw)


def bytes_to_tempfile(blob: bytes) -> str:
    """
    Decompress gzip-compressed and write it to temporary file.
    """
    pdb_text = gzip.decompress(blob).decode("utf-8")
    fd, path = tempfile.mkstemp(suffix=".pdb")
    with os.fdopen(fd, "w") as fh:
        fh.write(pdb_text)
    return path


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
 
 
def ddi_pair_to_bytes(pdb_file: str, chain_id_a: str, start_a: int, end_a: int, chain_id_b: str, start_b: int, end_b: int) -> bytes:
    """
    Slice out a single interacting domain pair (both domains of one DDI
    instance) from a full predicted complex PDB, keeping their original
    chain IDs.  Returns gzip-compressed PDB text bytes for storage in
    ddi_structure.pdb_gz.
 
    Unlike domain_to_bytes (single domain, single chain), this keeps two
    disjoint chain/residue windows in the same output structure so the
    pairwise contact geometry used in scoring (which residue pair on
    chain A contacts which on chain B) is preserved.
    """
    parser    = PDBParser(QUIET=True)
    structure = parser.get_structure("ddi_pair", pdb_file)
 
    buf    = io.StringIO()
    io_obj = PDBIO()
    io_obj.set_structure(structure)
    io_obj.save(
        buf,
        DdiPairSelect(chain_id_a, start_a, end_a, chain_id_b, start_b, end_b),
    )
 
    return gzip.compress(buf.getvalue().encode("utf-8"))



# ---------------------------------------------------------------------------
# Database helpers
# ---------------------------------------------------------------------------

def connect_db(db_path: str):
    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA synchronous  = OFF")
    conn.execute("PRAGMA journal_mode = MEMORY")
    return conn



# def add_pdb_to_mapping(conn, domain_id, protein_id, pdb_gz, source):
#     # Add single-domain psb to domain_protein_map table (for single-domain structures, not DDIs)
#     # If source is 'AF3', add to column pdb_gz_af, if source is 'RF2', add to column pdb_gz_rf, otherwise skip
#     if source == 'AF3':
#         conn.execute("""
#             UPDATE domain_protein_map
#             SET pdb_gz_af = ?
#             WHERE domain_id = ? AND protein_id = ?
#         """, (pdb_gz, domain_id, protein_id))
#     elif source == 'RF2':
#         conn.execute("""
#             UPDATE domain_protein_map
#             SET pdb_gz_rf = ?
#             WHERE domain_id = ? AND protein_id = ?
#         """, (pdb_gz, domain_id, protein_id))
#     else:
#         print(f"Warning: unknown source '{source}' for domain {domain_id}, protein {protein_id}. Skipping PDB addition.")


def get_ppis(conn):
    """Return all PPIs as (protein_id_a, uniprot_id_a, sequence_a,
                           protein_id_b, uniprot_id_b, sequence_b)."""
    return conn.execute("""
        SELECT
            p1.id         AS protein_id_a,
            p1.uniprot_id AS uniprot_id_a,
            p1.sequence   AS sequence_a,
            p2.id         AS protein_id_b,
            p2.uniprot_id AS uniprot_id_b,
            p2.sequence   AS sequence_b
        FROM protein_protein_interaction ppi
        JOIN protein p1 ON ppi.protein_id_a = p1.id
        JOIN protein p2 ON ppi.protein_id_b = p2.id
    """).fetchall()


def get_domain_mapping(conn, protein_id: int):
    """Return [(domain_id, start_pos, end_pos)] for a protein."""
    return conn.execute("""
        SELECT domain_id, start_pos, end_pos
        FROM domain_protein_map
        WHERE protein_id = ?
    """, (protein_id,)).fetchall()



def store_domain_slice(conn, ddi_id: int, domain_id_1: int, domain_id_2: int, protein_id_1: int, protein_id_2: int, pdb_gz: bytes, source_model: str) -> None:
    return conn.execute("""
        INSERT OR REPLACE INTO domain_structure
            (ddi_id, domain1, domain2, protein1, protein2, source, pdb_gz, z_score)
        VALUES (?, ?, ?, ?, ?, ?, ?, NULL)
    """, (ddi_id, domain_id_1, domain_id_2, protein_id_1, protein_id_2, source_model, pdb_gz))



def update_score(conn, ds_id: int, z_score: float) -> None:
    """Update the z_score for a domain_structure entry."""
    conn.execute("""
        UPDATE domain_structure
        SET z_score = ?
        WHERE id = ?
    """, (z_score, ds_id))
    conn.commit()



def check_ddi_exists(conn, domain_id_a: int, domain_id_b: int):
    """Return the ID of the DDI if it exists for this domain pair, otherwise None."""
    row = conn.execute("""
        SELECT id
        FROM domain_domain_interaction ddi
        WHERE ddi.domain_id_a = ? AND ddi.domain_id_b = ?
    """, (domain_id_a, domain_id_b)).fetchone()
    row_rev = conn.execute("""
        SELECT id
        FROM domain_domain_interaction ddi
        WHERE ddi.domain_id_a = ? AND ddi.domain_id_b = ?
    """, (domain_id_b, domain_id_a)).fetchone()

    id = row[0] if row else (row_rev[0] if row_rev else None)
    
    return id



#---------------------------------------------------------------------------
# Scoring Helper Functions
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

def residues_contact(resA, resB):
    """
    Return True if any atom pair between resA and resB satisfies the
    3did contact definition (H-bond, salt bridge, or vdW).
    """
    for atomA in resA.get_atoms():
        elemA = atomA.element.strip().upper() if atomA.element else ""
        for atomB in resB.get_atoms():
            elemB = atomB.element.strip().upper() if atomB.element else ""
            dist = atomA - atomB
            pair = tuple(sorted([elemA, elemB]))
            if pair == ('C', 'C') and dist <= CC_VDW:
                return True
            if pair == ('N', 'O') and dist <= NO_SB:
                return True
    return False



def get_domain_structures(db_path, source = None):
    """
    Get domain structure information from the SQLite database for all DDIs from 3DID.
    Returns a list of tuples: (ds_id, ddi_id, pdb_gz, source)
    """

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA synchronous  = OFF")
    conn.execute("PRAGMA journal_mode = MEMORY")


    ddis_positive = conn.execute("""
        SELECT id
        FROM domain_domain_interaction
        WHERE negative = 0
    """).fetchall()


    # Add chain information from complex_chain_map and pdb_gz from domain_structure
    domain_structures = []
    for ddi_id in ddis_positive:
        if not source:
        
            rows = conn.execute("""
                SELECT ds.id, ds.pdb_gz, ds.source
                FROM domain_structure ds
                WHERE ds.ddi_id = ?
            """, (ddi_id,)).fetchall()

        else:
            rows = conn.execute("""
                SELECT ds.id, ds.pdb_gz, ds.source
                FROM domain_structure ds
                WHERE ds.ddi_id = ? AND ds.source = ?
            """, (ddi_id[0], source)).fetchall()

        if not rows:
            continue


        for row in rows:
            ds_id, pdb_gz, source = row
            domain_structures.append((ds_id, ddi_id, pdb_gz, source))


    conn.close()
    return domain_structures



def get_domain_structures_for_scoring(db_path, source = None):
    """
    Get domain structure information from the SQLite database for all DDIs from 3DID.
    Returns a list of tuples: (ds_id, ddi_id, pdb_gz, source)
    """

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA synchronous  = OFF")


    domain_structures = conn.execute("""
        SELECT ds.id, ds.ddi_id, ds.pdb_gz, ds.source
        FROM domain_structure ds
        WHERE ds.source = ?
    """, (source,)).fetchall()


    conn.close()
    return domain_structures



#----------------------------------------------------------------------------
# Random PDB generation for AlphaFold
#----------------------------------------------------------------------------

# _BACKBONE_ATOMS = ["N", "CA", "C", "O"]
# _BACKBONE_OFFSETS = {  # offsets (A) from the residue's CA, idealised geometry
#     "N":  (-0.5, 1.4, 0.0),
#     "CA": (0.0, 0.0, 0.0),
#     "C":  (1.4, -0.5, 0.0),
#     "O":  (2.0, -1.5, 0.0),
# }
 
# _THREE_LETTER = {
#     'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
#     'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
#     'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
#     'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
# }
 

# def _write_mock_chain_atoms(lines: list, chain_id: str, sequence: str,
#                              atom_serial: list, x_offset: float) -> None:
#     """Append PDB ATOM lines for one chain with a simple backbone-only
#     random-walk trace, starting at x_offset along the X axis."""
#     rng = random.Random(f"{chain_id}-{sequence[:8]}-{len(sequence)}")
#     x, y, z = x_offset, 0.0, 0.0
#     rise = 3.8  # approx CA-CA spacing along the synthetic chain
 
#     for res_idx, aa in enumerate(sequence, start=1):
#         resname = _THREE_LETTER.get(aa.upper(), "GLY")
 
#         # small random perpendicular jitter so the structure isn't a
#         # perfectly straight (degenerate) line
#         jitter_y = rng.uniform(-1.0, 1.0)
#         jitter_z = rng.uniform(-1.0, 1.0)
#         cx, cy, cz = x, y + jitter_y, z + jitter_z
 
#         for atom_name in _BACKBONE_ATOMS:
#             dx, dy, dz = _BACKBONE_OFFSETS[atom_name]
#             ax, ay, az = cx + dx, cy + dy, cz + dz
#             element = atom_name[0]
#             atom_serial[0] += 1
#             lines.append(
#                 f"ATOM  {atom_serial[0]:5d}  {atom_name:<3s}{resname:>3s} "
#                 f"{chain_id}{res_idx:4d}    "
#                 f"{ax:8.3f}{ay:8.3f}{az:8.3f}{1.00:6.2f}{0.00:6.2f}"
#                 f"          {element:>2s}\n"
#             )
 
#         x += rise



# def mock_predict_complex(sequence_a: str, sequence_b: str, out_cif_dir: str,
#                           basename: str) -> str:
#     """
#     Generate a synthetic 'model_0' PDB structure for two sequences placed
#     on chains A and B, mimicking the file an AF3 run would eventually
#     produce (minus the CIF step -- this writes PDB directly).
 
#     This is NOT a real structure prediction. Coordinates are a randomised
#     backbone-only trace with no physical relevance -- it exists purely so
#     downstream code (domain slicing, contact geometry, scoring) can be
#     exercised before AF3 is wired up.
 
#     Returns the path to the written PDB file:
#         {out_cif_dir}/{basename}_model_0.pdb
#     """
#     pdb_path = os.path.join(out_cif_dir, f"{basename}_model_0.pdb")
 
#     lines = ["HEADER    MOCK AF3 PREDICTION (SYNTHETIC, FOR TESTING ONLY)\n"]
#     atom_serial = [0]
#     _write_mock_chain_atoms(lines, "A", sequence_a, atom_serial, x_offset=0.0)
#     lines.append("TER\n")
#     _write_mock_chain_atoms(lines, "B", sequence_b, atom_serial, x_offset=20.0)
#     lines.append("TER\n")
#     lines.append("END\n")
 
#     with open(pdb_path, "w") as fh:
#         fh.writelines(lines)
 
#     return pdb_path



# def mock_predict_complex_rf(sequence_a: str, sequence_b: str, out_cif_dir: str) -> str:
#     """
#     Generate a synthetic 'model_0' PDB structure for two sequences placed
#     on chains A and B, mimicking the file an AF3 run would eventually
#     produce (minus the CIF step -- this writes PDB directly).
 
#     This is NOT a real structure prediction. Coordinates are a randomised
#     backbone-only trace with no physical relevance -- it exists purely so
#     downstream code (domain slicing, contact geometry, scoring) can be
#     exercised before AF3 is wired up.
 
#     Returns the path to the written PDB file:
#         {out_cif_dir}/{basename}_model_0.pdb
#     """
#     pdb_path = os.path.join(out_cif_dir, f"model_0.pdb")
 
#     lines = ["HEADER    MOCK AF3 PREDICTION (SYNTHETIC, FOR TESTING ONLY)\n"]
#     atom_serial = [0]
#     _write_mock_chain_atoms(lines, "A", sequence_a, atom_serial, x_offset=0.0)
#     lines.append("TER\n")
#     _write_mock_chain_atoms(lines, "B", sequence_b, atom_serial, x_offset=20.0)
#     lines.append("TER\n")
#     lines.append("END\n")
 
#     with open(pdb_path, "w") as fh:
#         fh.writelines(lines)
 
#     return pdb_path