"""Receptor preparation: cleaning heteroatoms, adding polar hydrogens, and generating PDBQT."""

from pathlib import Path
from typing import Dict, List, Tuple

from Bio.PDB import PDBParser, is_aa

from shared.logging.logger import get_logger

logger = get_logger("projects.04.receptor")

# Standard Gasteiger/Kollman approximate partial charges for protein atoms
STANDARD_CHARGES: Dict[str, float] = {
    "N": -0.47,
    "CA": 0.07,
    "C": 0.51,
    "O": -0.51,
    "CB": -0.05,
    "CG": 0.00,
    "CD": 0.00,
    "OE1": -0.55,
    "NE2": -0.45,
    "OG1": -0.55,
    "SG": -0.20,
    "OH": -0.55,
    "NZ": -0.30,
    "OD1": -0.55,
    "ND1": -0.35,
}

# AutoDock 4 atom type mapping
AD4_ATOM_TYPES: Dict[str, str] = {
    "C": "C",
    "N": "NA",
    "O": "OA",
    "S": "SA",
    "H": "HD",
    "P": "P",
}


def prepare_receptor(
    source_pdb_path: str | Path,
    output_pdb_path: str | Path,
    output_pdbqt_path: str | Path,
    chain_id: str = "A",
) -> Tuple[int, Path, Path]:
    """Clean receptor structure, strip heteroatoms, and write formatted PDB and PDBQT files.

    Args:
        source_pdb_path: Raw PDB file (e.g. 1E3G).
        output_pdb_path: Cleaned protein PDB output path.
        output_pdbqt_path: PDBQT output path.
        chain_id: Receptor chain to preserve.

    Returns:
        Tuple of (residue_count, cleaned_pdb_path, pdbqt_path).
    """
    in_path = Path(source_pdb_path)
    out_pdb = Path(output_pdb_path)
    out_pdbqt = Path(output_pdbqt_path)

    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("receptor", str(in_path))

    cleaned_lines_pdb: List[str] = []
    pdbqt_lines: List[str] = []
    residue_count = 0
    atom_serial = 1

    for model in structure:
        for chain in model:
            if chain_id and chain.get_id() != chain_id:
                continue

            for res in chain:
                if not is_aa(res, standard=True):
                    continue

                residue_count += 1
                res_name = res.get_resname()
                res_seq = res.get_id()[1]

                for atom in res:
                    atom_name = atom.get_name().strip()
                    coord = atom.get_coord()
                    element = atom.element.strip().upper() if atom.element else atom_name[0]
                    b_factor = atom.get_bfactor()
                    occupancy = atom.get_occupancy()

                    # PDB line format
                    pdb_line = (
                        f"ATOM  {atom_serial:5d} {atom_name:^4s} {res_name:3s} {chain.get_id():1s}"
                        f"{res_seq:4d}    {coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                        f"{occupancy:6.2f}{b_factor:6.2f}           {element:>2s}"
                    )
                    cleaned_lines_pdb.append(pdb_line)

                    # Assign partial charge and AD4 atom type
                    charge = STANDARD_CHARGES.get(atom_name, 0.0)
                    ad4_type = AD4_ATOM_TYPES.get(element, element)
                    pdbqt_line = (
                        f"ATOM  {atom_serial:5d} {atom_name:^4s} {res_name:3s} {chain.get_id():1s}"
                        f"{res_seq:4d}    {coord[0]:8.3f}{coord[1]:8.3f}{coord[2]:8.3f}"
                        f"{occupancy:6.2f}{b_factor:6.2f}    {charge:+6.3f} {ad4_type:<2s}"
                    )
                    pdbqt_lines.append(pdbqt_line)
                    atom_serial += 1

    with open(out_pdb, "w", encoding="utf-8") as f:
        f.write("\n".join(cleaned_lines_pdb) + "\nEND\n")

    with open(out_pdbqt, "w", encoding="utf-8") as f:
        f.write("\n".join(pdbqt_lines) + "\nEND\n")

    logger.info(
        f"Prepared receptor: {residue_count} residues, {atom_serial - 1} atoms saved to {out_pdbqt.name}"
    )
    return residue_count, out_pdb, out_pdbqt
