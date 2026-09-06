"""Protein-ligand interaction profiling (hydrogen bonds, hydrophobic contacts)."""

from pathlib import Path
from typing import Any, Dict, List, Tuple

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, is_aa

from shared.logging.logger import get_logger

logger = get_logger("projects.04.interactions")


def parse_pose_atoms(
    pdbqt_file: str | Path, target_mode: int = 1
) -> List[Tuple[str, np.ndarray, str]]:
    """Parse atom names, 3D coordinates, and elements for a specific mode from a multi-model PDBQT."""
    path = Path(pdbqt_file)
    atoms: List[Tuple[str, np.ndarray, str]] = []
    in_mode = False

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(f"MODEL {target_mode}"):
                in_mode = True
            elif line.startswith("ENDMDL") and in_mode:
                break
            elif in_mode and (line.startswith("ATOM") or line.startswith("HETATM")):
                try:
                    name = line[12:16].strip()
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    element = line[76:78].strip().upper() if len(line) >= 78 else name[0]
                    atoms.append((name, np.array([x, y, z]), element))
                except Exception:
                    continue

    return atoms


def analyze_pose_interactions(
    receptor_pdb: str | Path,
    pose_pdbqt: str | Path,
    ligand_id: str,
    target_mode: int = 1,
    hbond_cutoff: float = 3.5,
    hydrophobic_cutoff: float = 4.0,
) -> pd.DataFrame:
    """Detect hydrogen bonds and hydrophobic contacts between a docked pose and the receptor."""
    parser = PDBParser(QUIET=True)
    receptor = parser.get_structure("rec", str(receptor_pdb))

    lig_atoms = parse_pose_atoms(pose_pdbqt, target_mode=target_mode)
    if not lig_atoms:
        logger.warning(f"No atoms parsed for {ligand_id} mode {target_mode}.")
        return pd.DataFrame()

    interactions: List[Dict[str, Any]] = []

    for model in receptor:
        for chain in model:
            for res in chain:
                if not is_aa(res, standard=True):
                    continue

                res_num = res.get_id()[1]
                res_name = res.get_resname()

                for patom in res:
                    pcoord = patom.get_coord()
                    pelem = patom.element.strip().upper() if patom.element else patom.get_name()[0]

                    for latom_name, lcoord, lelem in lig_atoms:
                        dist = float(np.linalg.norm(pcoord - lcoord))

                        # Hydrogen bond: O or N pairs within hbond_cutoff
                        if dist <= hbond_cutoff and pelem in ("O", "N") and lelem in ("O", "N"):
                            interactions.append(
                                {
                                    "ligand_id": ligand_id,
                                    "mode": target_mode,
                                    "residue_number": res_num,
                                    "residue_name": res_name,
                                    "interaction_type": "Hydrogen Bond",
                                    "distance_angstroms": round(dist, 2),
                                    "receptor_atom": patom.get_name().strip(),
                                    "ligand_atom": latom_name,
                                }
                            )

                        # Hydrophobic contact: C-C pairs within hydrophobic_cutoff
                        elif dist <= hydrophobic_cutoff and pelem == "C" and lelem in ("C", "A"):
                            interactions.append(
                                {
                                    "ligand_id": ligand_id,
                                    "mode": target_mode,
                                    "residue_number": res_num,
                                    "residue_name": res_name,
                                    "interaction_type": "Hydrophobic Contact",
                                    "distance_angstroms": round(dist, 2),
                                    "receptor_atom": patom.get_name().strip(),
                                    "ligand_atom": latom_name,
                                }
                            )

    df = pd.DataFrame(interactions)
    if not df.empty:
        # Keep closest contact per residue & interaction type
        df = df.sort_values("distance_angstroms").drop_duplicates(
            subset=["residue_number", "interaction_type"]
        )

    logger.info(f"Interaction analysis for {ligand_id}: {len(df)} contacts detected.")
    return df
