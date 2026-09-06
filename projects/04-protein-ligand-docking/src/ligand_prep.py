"""Ligand preparation and 3D conformer optimization with RDKit."""

from dataclasses import asdict, dataclass
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from shared.logging.logger import get_logger

logger = get_logger("projects.04.ligand")


@dataclass
class PreparedLigand:
    """Standardized representation of a prepared 3D ligand."""

    ligand_id: str
    name: str
    smiles: str
    formula: str
    molecular_weight: float
    num_rotatable_bonds: int
    h_bond_donors: int
    h_bond_acceptors: int
    logp: float
    pdbqt_path: str
    mol_block: str

    def to_dict(self) -> dict:
        return asdict(self)


def prepare_ligand_from_smiles(
    smiles: str,
    ligand_id: str,
    output_dir: str | Path,
    name: str = "",
    random_seed: int = 42,
) -> PreparedLigand:
    """Generate 3D conformer, add hydrogens, energy minimize with MMFF94, and export PDBQT.

    Args:
        smiles: Canonical or isomeric SMILES string.
        ligand_id: Short unique identifier (e.g. 'DHT').
        output_dir: Directory to save generated structure files.
        name: Human-readable name.
        random_seed: Random seed for deterministic conformer generation.

    Returns:
        PreparedLigand dataclass.
    """
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string for ligand '{ligand_id}': {smiles}")

    mol = Chem.AddHs(mol)

    # 3D Conformer Generation with ETKDGv3
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    res = AllChem.EmbedMolecule(mol, params)
    if res != 0:
        # Fallback to standard embedding
        AllChem.EmbedMolecule(mol, randomSeed=random_seed)

    # Energy Minimization using MMFF94 force field
    try:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception as e:
        logger.warning(f"MMFF optimization fallback for {ligand_id}: {e}")

    # Compute Gasteiger partial charges
    AllChem.ComputeGasteigerCharges(mol)

    # Compute Cheminformatics Descriptors
    mw = round(Descriptors.MolWt(mol), 2)
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
    rot_bonds = Descriptors.NumRotatableBonds(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    logp = round(Descriptors.MolLogP(mol), 2)

    # Write PDBQT format
    pdbqt_path = out_dir / f"{ligand_id}.pdbqt"
    conf = mol.GetConformer()
    pdbqt_lines = [
        f"REMARK  LIGAND {ligand_id} ({name})",
        f"REMARK  FORMULA: {formula}  MW: {mw}  ROTATABLE_BONDS: {rot_bonds}",
        "ROOT",
    ]

    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        elem = atom.GetSymbol().upper()
        charge = (
            float(atom.GetProp("_GasteigerCharge")) if atom.HasProp("_GasteigerCharge") else 0.0
        )
        # AD4 atom type
        ad4_type = elem
        if elem == "C" and atom.GetIsAromatic():
            ad4_type = "A"
        elif elem == "N":
            ad4_type = "NA" if hba > 0 else "N"
        elif elem == "O":
            ad4_type = "OA"

        line = (
            f"ATOM  {atom.GetIdx() + 1:5d} {elem:^4s} UNK A   1    "
            f"{pos.x:8.3f}{pos.y:8.3f}{pos.z:8.3f}  1.00  0.00    "
            f"{charge:+6.3f} {ad4_type:<2s}"
        )
        pdbqt_lines.append(line)

    pdbqt_lines.append("ENDROOT")
    pdbqt_lines.append(f"TORSDOF {rot_bonds}")

    with open(pdbqt_path, "w", encoding="utf-8") as f:
        f.write("\n".join(pdbqt_lines) + "\n")

    mol_block = Chem.MolToMolBlock(mol)
    logger.info(f"Prepared ligand {ligand_id}: MW={mw}, Formula={formula}, PDBQT={pdbqt_path.name}")

    return PreparedLigand(
        ligand_id=ligand_id,
        name=name or ligand_id,
        smiles=smiles,
        formula=formula,
        molecular_weight=mw,
        num_rotatable_bonds=rot_bonds,
        h_bond_donors=hbd,
        h_bond_acceptors=hba,
        logp=logp,
        pdbqt_path=str(pdbqt_path),
        mol_block=mol_block,
    )
