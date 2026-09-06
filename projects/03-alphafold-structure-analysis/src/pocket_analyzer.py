"""Ligand-binding pocket identification and active-site residue evaluation."""

from pathlib import Path
from typing import Any, Dict, List, Set

import numpy as np
import pandas as pd
from Bio.PDB import NeighborSearch, PDBParser, is_aa

from shared.logging.logger import get_logger

logger = get_logger("projects.03.pocket")


def find_pocket_residues(
    pdb_path: str | Path,
    ligand_resname: str = "DHT",
    distance_cutoff: float = 4.5,
) -> tuple[pd.DataFrame, Set[int]]:
    """Identify all protein amino acid residues with an atom within distance_cutoff of the specified ligand.

    Args:
        pdb_path: Path to experimental complex PDB (e.g. 1E3G).
        ligand_resname: 3-letter code of ligand (e.g. 'DHT').
        distance_cutoff: Spatial contact threshold in Angstroms.

    Returns:
        Tuple of (pocket_residues_df, set_of_residue_numbers).
    """
    path = Path(pdb_path)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", str(path))

    # Collect ligand atoms
    ligand_atoms = []
    for model in structure:
        for chain in model:
            for res in chain:
                if res.get_resname() == ligand_resname:
                    ligand_atoms.extend(res.get_atoms())

    if not ligand_atoms:
        logger.warning(f"No ligand residues matching '{ligand_resname}' found in {path.name}.")
        return pd.DataFrame(), set()

    # Collect all protein atoms
    protein_atoms = []
    for model in structure:
        for chain in model:
            for res in chain:
                if is_aa(res, standard=True):
                    protein_atoms.extend(res.get_atoms())

    ns = NeighborSearch(protein_atoms)
    contact_residues = set()
    records: List[Dict[str, Any]] = []

    for latom in ligand_atoms:
        nearby_atoms = ns.search(latom.get_coord(), distance_cutoff, level="A")
        for patom in nearby_atoms:
            parent_res = patom.get_parent()
            res_id = parent_res.get_id()[1]
            dist = float(np.linalg.norm(latom.get_coord() - patom.get_coord()))
            contact_residues.add(res_id)
            records.append(
                {
                    "residue_number": res_id,
                    "residue_name": parent_res.get_resname(),
                    "chain": parent_res.get_parent().get_id(),
                    "atom_name": patom.get_name(),
                    "ligand_atom": latom.get_name(),
                    "distance_angstroms": round(dist, 2),
                }
            )

    df = pd.DataFrame(records)
    # Deduplicate to find minimum distance per residue
    if not df.empty:
        summary_df = (
            df.groupby(["residue_number", "residue_name", "chain"])
            .agg({"distance_angstroms": "min"})
            .reset_index()
            .sort_values("distance_angstroms")
        )
    else:
        summary_df = pd.DataFrame(
            columns=["residue_number", "residue_name", "chain", "distance_angstroms"]
        )

    logger.info(
        f"Identified {len(summary_df)} binding pocket residues within {distance_cutoff} Å of {ligand_resname}."
    )
    return summary_df, contact_residues
