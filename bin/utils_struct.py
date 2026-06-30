#! /usr/bin/env python3

import gzip
import io
import os
import tempfile
import sqlite3
import random

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.PDBIO import PDBIO, Select


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


def create_tables(conn):
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS complex_chain_map (
            id           INTEGER PRIMARY KEY AUTOINCREMENT,
            protein_id_a INTEGER NOT NULL REFERENCES protein(id),
            protein_id_b INTEGER NOT NULL REFERENCES protein(id),
            source       TEXT    NOT NULL, -- 'AF3' | 'RF2'
            chain_id_a   TEXT    NOT NULL,
            chain_id_b   TEXT    NOT NULL,
            UNIQUE (protein_id_a, protein_id_b, source)
        );

        CREATE TABLE IF NOT EXISTS domain_structure (
            id         INTEGER PRIMARY KEY AUTOINCREMENT,
            ddi_id INTEGER NOT NULL REFERENCES domain_domain_interaction(id),
            protein1 INTEGER NOT NULL REFERENCES protein(id),
            protein2 INTEGER NOT NULL REFERENCES protein(id),
            source     TEXT    NOT NULL,
            pdb_gz     BLOB    NOT NULL,
            z_score    REAL,
            UNIQUE (ddi_id, protein1, protein2, source)
        );
    """)
    conn.commit()


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


def store_chain_map(conn, protein_id_a: int, protein_id_b: int, chain_id_a: str, chain_id_b: str, source_model: str) -> None:
    return conn.execute("""
        INSERT OR REPLACE INTO complex_chain_map
            (protein_id_a, protein_id_b, source, chain_id_a, chain_id_b)
        VALUES (?, ?, ?, ?, ?)
    """, (protein_id_a, protein_id_b, source_model, chain_id_a, chain_id_b))


def store_domain_slice(conn, protein_id_1: int, protein_id_2: int, ddi_id: int, pdb_gz: bytes, source_model: str) -> None:
    return conn.execute("""
        INSERT OR REPLACE INTO domain_structure
            (ddi_id, protein1, protein2,  source, pdb_gz, z_score)
        VALUES (?, ?, ?, ?, ?, NULL)
    """, (ddi_id, protein_id_1, protein_id_2, source_model, pdb_gz))


def check_ddi_exists(conn, domain_id_a: int, domain_id_b: int):
    """Return the ID of the DDI if it exists for this domain pair, otherwise None."""
    row = conn.execute("""
        SELECT id
        FROM domain_domain_interaction ddi
        WHERE ddi.domain1 = ? AND ddi.domain2 = ?
    """, (domain_id_a, domain_id_b)).fetchone()
    row_rev = conn.execute("""
        SELECT id
        FROM domain_domain_interaction ddi
        WHERE ddi.domain1 = ? AND ddi.domain2 = ?
    """, (domain_id_b, domain_id_a)).fetchone()

    id = row[0] if row else (row_rev[0] if row_rev else None)
    
    return id



def get_predicted_complexes(conn, source: str):
    return conn.execute("""
        SELECT
            p1.id AS protein_id_a, p1.uniprot_id AS uniprot_id_a, ccm.chain_id_a,
            p2.id AS protein_id_b, p2.uniprot_id AS uniprot_id_b, ccm.chain_id_b
        FROM complex_chain_map ccm
        JOIN protein p1 ON ccm.protein_id_a = p1.id
        JOIN protein p2 ON ccm.protein_id_b = p2.id
        WHERE ccm.source = ?
    """, (source,)).fetchall()



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



def get_domain_structures(db_path):
    """
    Get domain structure information from the SQLite database for all DDIs from 3DID.
    Returns a list of tuples: (ddi_id, protein_id_a, protein_id_b, chain_id_a, chain_id_b, pdb_gz)
    """

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA foreign_keys = ON")
    conn.execute("PRAGMA synchronous  = OFF")
    conn.execute("PRAGMA journal_mode = MEMORY")


    ddis_from_3did = conn.execute("""
        SELECT id, domain_id_a, domain_id_b
        FROM domain_domain_interaction
        WHERE source = '3DID'
    """).fetchall()

    # Add chain information from complex_chain_map and pdb_gz from domain_structure
    domain_structures = []
    for ddi_id, domain_id_a, domain_id_b in ddis_from_3did:
        row = conn.execute("""
            SELECT d2.id, ds.protein1, ds.protein2, ccm.chain_id_a, ccm.chain_id_b, ds.pdb_gz
            FROM domain_structure ds
            JOIN complex_chain_map ccm ON ccm.protein_id_a = ds.protein1 AND ccm.protein_id_b = ds.protein2 AND ccm.source = ds.source
            WHERE ds.ddi_id = ?
        """, (ddi_id,)).fetchone()
        
        if row:
            ds_id, protein_id_a, protein_id_b, chain_id_a, chain_id_b, pdb_gz = row
            domain_structures.append((ds_id, ddi_id, chain_id_a, chain_id_b, pdb_gz))


    conn.close()
    return domain_structures



#----------------------------------------------------------------------------
# Random PDB generation for AlphaFold
#----------------------------------------------------------------------------

_BACKBONE_ATOMS = ["N", "CA", "C", "O"]
_BACKBONE_OFFSETS = {  # offsets (A) from the residue's CA, idealised geometry
    "N":  (-0.5, 1.4, 0.0),
    "CA": (0.0, 0.0, 0.0),
    "C":  (1.4, -0.5, 0.0),
    "O":  (2.0, -1.5, 0.0),
}
 
_THREE_LETTER = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL',
}
 

def _write_mock_chain_atoms(lines: list, chain_id: str, sequence: str,
                             atom_serial: list, x_offset: float) -> None:
    """Append PDB ATOM lines for one chain with a simple backbone-only
    random-walk trace, starting at x_offset along the X axis."""
    rng = random.Random(f"{chain_id}-{sequence[:8]}-{len(sequence)}")
    x, y, z = x_offset, 0.0, 0.0
    rise = 3.8  # approx CA-CA spacing along the synthetic chain
 
    for res_idx, aa in enumerate(sequence, start=1):
        resname = _THREE_LETTER.get(aa.upper(), "GLY")
 
        # small random perpendicular jitter so the structure isn't a
        # perfectly straight (degenerate) line
        jitter_y = rng.uniform(-1.0, 1.0)
        jitter_z = rng.uniform(-1.0, 1.0)
        cx, cy, cz = x, y + jitter_y, z + jitter_z
 
        for atom_name in _BACKBONE_ATOMS:
            dx, dy, dz = _BACKBONE_OFFSETS[atom_name]
            ax, ay, az = cx + dx, cy + dy, cz + dz
            element = atom_name[0]
            atom_serial[0] += 1
            lines.append(
                f"ATOM  {atom_serial[0]:5d}  {atom_name:<3s}{resname:>3s} "
                f"{chain_id}{res_idx:4d}    "
                f"{ax:8.3f}{ay:8.3f}{az:8.3f}{1.00:6.2f}{0.00:6.2f}"
                f"          {element:>2s}\n"
            )
 
        x += rise



def mock_predict_complex(sequence_a: str, sequence_b: str, out_cif_dir: str,
                          basename: str) -> str:
    """
    Generate a synthetic 'model_0' PDB structure for two sequences placed
    on chains A and B, mimicking the file an AF3 run would eventually
    produce (minus the CIF step -- this writes PDB directly).
 
    This is NOT a real structure prediction. Coordinates are a randomised
    backbone-only trace with no physical relevance -- it exists purely so
    downstream code (domain slicing, contact geometry, scoring) can be
    exercised before AF3 is wired up.
 
    Returns the path to the written PDB file:
        {out_cif_dir}/{basename}_model_0.pdb
    """
    pdb_path = os.path.join(out_cif_dir, f"{basename}_model_0.pdb")
 
    lines = ["HEADER    MOCK AF3 PREDICTION (SYNTHETIC, FOR TESTING ONLY)\n"]
    atom_serial = [0]
    _write_mock_chain_atoms(lines, "A", sequence_a, atom_serial, x_offset=0.0)
    lines.append("TER\n")
    _write_mock_chain_atoms(lines, "B", sequence_b, atom_serial, x_offset=20.0)
    lines.append("TER\n")
    lines.append("END\n")
 
    with open(pdb_path, "w") as fh:
        fh.writelines(lines)
 
    return pdb_path



