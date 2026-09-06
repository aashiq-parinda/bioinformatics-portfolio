"""PDB and structural format parser and atom extractor."""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, List

import numpy as np
from Bio.PDB import PDBParser, Structure, is_aa
from Bio.PDB.PDBIO import PDBIO, Select

from shared.logging.logger import get_logger

logger = get_logger("shared.io.pdb")


@dataclass
class PDBResidue:
    """Residue information with sequence position and B-factor / pLDDT."""

    res_id: int
    res_name: str
    chain_id: str
    b_factor_mean: float
    is_standard_aa: bool


class NotHeteroSelect(Select):
    """Filter class to retain only standard protein amino acid residues."""

    def accept_residue(self, residue: Any) -> int:
        return 1 if is_aa(residue, standard=True) else 0


def parse_pdb_structure(pdb_path: str | Path, structure_id: str = "prot") -> Structure.Structure:
    """Parse a PDB file using Biopython PDBParser."""
    path = Path(pdb_path)
    if not path.is_file():
        raise FileNotFoundError(f"PDB file not found: {path}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(structure_id, str(path))
    logger.info(f"Loaded structure '{structure_id}' from {path.name}")
    return structure


def extract_residues(structure: Structure.Structure) -> List[PDBResidue]:
    """Extract ordered residue summary including mean B-factors/pLDDT from a structure."""
    residues: List[PDBResidue] = []
    for model in structure:
        for chain in model:
            for res in chain:
                b_factors = [atom.get_bfactor() for atom in res.get_atoms()]
                mean_b = float(np.mean(b_factors)) if b_factors else 0.0
                residues.append(
                    PDBResidue(
                        res_id=res.get_id()[1],
                        res_name=res.get_resname(),
                        chain_id=chain.get_id(),
                        b_factor_mean=round(mean_b, 2),
                        is_standard_aa=is_aa(res, standard=True),
                    )
                )
    return residues


def clean_pdb(input_path: str | Path, output_path: str | Path) -> None:
    """Extract only standard protein amino acid residues, strip heteroatoms and water molecules."""
    structure = parse_pdb_structure(input_path)
    io = PDBIO()
    io.set_structure(structure)
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    io.save(str(out), select=NotHeteroSelect())
    logger.info(f"Cleaned protein saved to {out}")
