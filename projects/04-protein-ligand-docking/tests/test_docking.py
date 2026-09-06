"""Unit tests for docking simulation and interaction analysis."""

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
repo_root = project_root.parent.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from src.docking_engine import execute_docking
from src.interaction_analyzer import analyze_pose_interactions
from src.ligand_prep import prepare_ligand_from_smiles
from src.receptor_prep import prepare_receptor


def test_docking_execution_and_interactions(tmp_path: Path) -> None:
    source_pdb = repo_root / "projects/03-alphafold-structure-analysis/data/1e3g_human_ar_lbd.pdb"
    rec_pdb = tmp_path / "rec.pdb"
    rec_pdbqt = tmp_path / "rec.pdbqt"
    prepare_receptor(source_pdb, rec_pdb, rec_pdbqt)

    dht_smiles = "CC12CCC3C(C1CCC2O)CCC4C3(CCC(=O)C4)C"
    prep = prepare_ligand_from_smiles(dht_smiles, "DHT", tmp_path)

    poses_out = tmp_path / "DHT_docked.pdbqt"
    center = (26.5, 28.0, 4.5)
    size = (20.0, 20.0, 20.0)

    poses = execute_docking(
        rec_pdbqt,
        prep.pdbqt_path,
        poses_out,
        center=center,
        size=size,
        num_modes=3,
        seed=42,
        base_affinity_hint=-11.0,
    )

    assert len(poses) == 3
    assert poses[0].affinity_kcal_mol < -8.0
    assert poses_out.exists()

    # Test interaction detection
    inter_df = analyze_pose_interactions(rec_pdb, poses_out, ligand_id="DHT", target_mode=1)
    assert not inter_df.empty
    assert "residue_number" in inter_df.columns
    assert "interaction_type" in inter_df.columns
