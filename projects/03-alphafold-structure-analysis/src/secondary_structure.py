"""Secondary structure classification via backbone dihedral angles."""

import math
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
from Bio.PDB import PDBParser, Polypeptide, is_aa

from shared.logging.logger import get_logger

logger = get_logger("projects.03.secstruct")


def compute_secondary_structure(pdb_path: str | Path) -> tuple[pd.DataFrame, Dict[str, float]]:
    """Compute per-residue secondary structure assignments using backbone phi/psi dihedral angles."""
    path = Path(pdb_path)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(path))

    records: List[Dict[str, Any]] = []

    for model in structure:
        for chain in model:
            polypeptides = Polypeptide.PPBuilder().build_peptides(chain)
            for pp in polypeptides:
                phi_psi_list = pp.get_phi_psi_list()
                for res, (phi, psi) in zip(pp, phi_psi_list, strict=False):
                    if not is_aa(res, standard=True):
                        continue

                    phi_deg = math.degrees(phi) if phi is not None else None
                    psi_deg = math.degrees(psi) if psi is not None else None

                    # Ramachandran region classification
                    ss_type = "Loop / Coil"
                    if phi_deg is not None and psi_deg is not None:
                        if -160 <= phi_deg <= -35 and -70 <= psi_deg <= 50:
                            ss_type = "Alpha-Helix"
                        elif -180 <= phi_deg <= -50 and (80 <= psi_deg <= 180 or psi_deg <= -170):
                            ss_type = "Beta-Sheet"
                        elif 40 <= phi_deg <= 100 and -40 <= psi_deg <= 60:
                            ss_type = "Left-handed Helix"

                    records.append(
                        {
                            "residue_number": res.get_id()[1],
                            "residue_name": res.get_resname(),
                            "chain": chain.get_id(),
                            "phi": round(phi_deg, 1) if phi_deg else None,
                            "psi": round(psi_deg, 1) if psi_deg else None,
                            "secondary_structure": ss_type,
                        }
                    )

    df = pd.DataFrame(records)
    total = len(df)
    if total == 0:
        return df, {"Alpha-Helix": 0.0, "Beta-Sheet": 0.0, "Loop / Coil": 100.0}

    counts = df["secondary_structure"].value_counts()
    percentages = {
        ss: round(float(counts.get(ss, 0) / total * 100.0), 1)
        for ss in ["Alpha-Helix", "Beta-Sheet", "Loop / Coil"]
    }

    logger.info(
        f"Secondary structure: {percentages['Alpha-Helix']}% Helix, "
        f"{percentages['Beta-Sheet']}% Sheet, {percentages['Loop / Coil']}% Coil"
    )
    return df, percentages
