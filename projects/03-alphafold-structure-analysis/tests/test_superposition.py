"""Unit tests for C-alpha superposition and RMSD calculation."""

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
repo_root = project_root.parent.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from src.structure_comparison import calculate_ca_rmsd


def test_calculate_ca_rmsd() -> None:
    ref_pdb = project_root / "data" / "1e3g_human_ar_lbd.pdb"
    target_pdb = project_root / "data" / "af_p10275_human_ar.pdb"

    # Superimpose AR LBD (residues 676-919)
    rmsd, n_atoms = calculate_ca_rmsd(ref_pdb, target_pdb, res_start=676, res_end=919)

    assert n_atoms > 150
    # AlphaFold LBD should agree well with the crystal structure (RMSD typically < 3.5 Å)
    assert 0.1 < rmsd < 5.0
