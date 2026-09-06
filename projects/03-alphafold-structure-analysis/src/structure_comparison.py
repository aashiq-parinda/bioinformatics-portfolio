"""Structural superposition and C-alpha RMSD computation."""

from pathlib import Path

from Bio.PDB import PDBParser, Superimposer, is_aa

from shared.logging.logger import get_logger

logger = get_logger("projects.03.superimpose")


def calculate_ca_rmsd(
    ref_pdb_path: str | Path,
    target_pdb_path: str | Path,
    res_start: int = 676,
    res_end: int = 919,
) -> tuple[float, int]:
    """Superimpose target PDB onto reference PDB across a residue range and compute C-alpha RMSD.

    Args:
        ref_pdb_path: Experimental reference structure (e.g. 1E3G).
        target_pdb_path: Predicted model structure (e.g. AlphaFold).
        res_start: Starting residue number.
        res_end: Ending residue number.

    Returns:
        Tuple of (RMSD_in_Angstroms, number_of_aligned_atoms).
    """
    parser = PDBParser(QUIET=True)
    ref_struct = parser.get_structure("ref", str(ref_pdb_path))
    target_struct = parser.get_structure("target", str(target_pdb_path))

    # Map C-alpha atoms by residue number
    ref_ca_dict = {}
    for res in ref_struct.get_residues():
        if is_aa(res, standard=True) and "CA" in res:
            res_num = res.get_id()[1]
            if res_start <= res_num <= res_end:
                ref_ca_dict[res_num] = res["CA"]

    target_ca_dict = {}
    for res in target_struct.get_residues():
        if is_aa(res, standard=True) and "CA" in res:
            res_num = res.get_id()[1]
            if res_start <= res_num <= res_end:
                target_ca_dict[res_num] = res["CA"]

    # Intersect common residues
    common_res = sorted(set(ref_ca_dict.keys()) & set(target_ca_dict.keys()))
    if len(common_res) < 10:
        logger.warning(
            f"Only {len(common_res)} common residues between structures. Returning 0.0 RMSD."
        )
        return 0.0, len(common_res)

    ref_atoms = [ref_ca_dict[r] for r in common_res]
    target_atoms = [target_ca_dict[r] for r in common_res]

    superimposer = Superimposer()
    superimposer.set_atoms(ref_atoms, target_atoms)
    rmsd = float(superimposer.rms)

    logger.info(
        f"C-alpha superposition over {len(common_res)} residues yielded RMSD: {rmsd:.2f} Å."
    )
    return round(rmsd, 2), len(common_res)
