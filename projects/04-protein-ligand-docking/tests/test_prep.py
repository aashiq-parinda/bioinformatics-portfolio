"""Unit tests for receptor and ligand preparation."""

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
repo_root = project_root.parent.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from src.ligand_prep import prepare_ligand_from_smiles
from src.receptor_prep import prepare_receptor


def test_prepare_ligand_dht(tmp_path: Path) -> None:
    # 5alpha-Dihydrotestosterone SMILES
    dht_smiles = "CC12CCC3C(C1CCC2O)CCC4C3(CCC(=O)C4)C"
    prep = prepare_ligand_from_smiles(dht_smiles, "DHT", tmp_path, name="Dihydrotestosterone")

    assert prep.ligand_id == "DHT"
    assert 285 < prep.molecular_weight < 295
    assert prep.formula == "C19H30O2"
    assert Path(prep.pdbqt_path).exists()
    assert Path(prep.pdbqt_path).stat().st_size > 500


def test_prepare_receptor(tmp_path: Path) -> None:
    source_pdb = repo_root / "projects/03-alphafold-structure-analysis/data/1e3g_human_ar_lbd.pdb"
    out_pdb = tmp_path / "target.pdb"
    out_pdbqt = tmp_path / "target.pdbqt"

    res_count, p_pdb, p_pdbqt = prepare_receptor(source_pdb, out_pdb, out_pdbqt, chain_id="A")

    assert res_count > 200
    assert out_pdb.exists()
    assert out_pdbqt.exists()
