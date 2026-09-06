"""Unit tests for active-site pocket residue identification."""

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
repo_root = project_root.parent.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from src.pocket_analyzer import find_pocket_residues


def test_find_pocket_residues_1e3g() -> None:
    crystal_pdb = project_root / "data" / "1e3g_human_ar_lbd.pdb"
    df, contact_set = find_pocket_residues(crystal_pdb, ligand_resname="R18", distance_cutoff=4.5)

    assert not df.empty
    assert len(contact_set) > 5

    # Known essential AR pocket residues should be detected
    # 705 (Asn), 877 (Thr), 764 (Phe), or 745 (Met)
    assert 705 in contact_set or 877 in contact_set or 764 in contact_set
