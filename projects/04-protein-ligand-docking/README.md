# Project 04: Reproducible Protein-Ligand Docking Pipeline

[![Workflow: Snakemake](https://img.shields.io/badge/Workflow-Snakemake-blue.svg)](workflows/Snakefile)
[![Docker: Supported](https://img.shields.io/badge/Docker-Ready-2496ED.svg)](workflows/Dockerfile)
[![Tests: Pytest](https://img.shields.io/badge/tests-passing-brightgreen.svg)](tests/)

A research-oriented, parameterized molecular docking workflow assessing non-covalent binding interactions, pose ranking, and contact fingerprinting of steroid receptor modulators complexed with the Human Androgen Receptor Ligand-Binding Domain (LBD).

---

## 🧬 Scientific Overview & Chemical Systems

The pipeline benchmarks three distinct small-molecule pharmacology classes against the Human AR-LBD pocket:
1. **$5\alpha$-Dihydrotestosterone (DHT)**: High-affinity natural endogenous agonist (`PubChem CID: 10635`).
2. **Enzalutamide (MDV3100)**: Non-steroidal 2nd-generation antagonist approved for castration-resistant prostate cancer (`PubChem CID: 15951529`).
3. **Testosterone**: Primary circulating androgen (`PubChem CID: 6010`).

---

## ⚠️ MANDATORY SCIENTIFIC & SAFETY DISCLAIMER

> **COMPUTATIONAL RESEARCH NOTICE:**
> Molecular docking is a computational hypothesis-generation and virtual screening methodology. Docking scores ($\Delta G$ in kcal/mol) **DO NOT** prove biological efficacy, clinical potency, therapeutic safety, or bioavailability. Under no circumstances should docking predictions be interpreted as clinical dosing protocols, medical advice, or bodybuilding recommendations. All findings require rigorous experimental wet-lab assay validation.

---

## 🔬 Computational Methodology

1. **Target Receptor Preparation**: Strip crystallographic water molecules and co-crystallized ligand from PDB `1E3G`, parameterize with Gasteiger partial charges, and format into PDBQT.
2. **Ligand Preparation with RDKit**: Ingest canonical SMILES, generate 3D conformers via ETKDGv3, energy minimize with the MMFF94 force field, compute partial charges, and export PDBQT.
3. **Docking Execution**: Parameterized search space centered on the binding cleft ($(26.5, 28.0, 4.5)$, $22\times 22\times 22$ Å box), random seed = 42, exhaustiveness = 8.
4. **Interaction Fingerprinting**: Automated detection of hydrogen bonds ($\le 3.5$ Å) and hydrophobic contact shells ($\le 4.0$ Å) against key catalytic residues (Asn705, Thr877, Phe764, Met745, Leu704).

---

## 🛠️ Execution Instructions

### Run via Python Pipeline
```bash
python3 projects/04-protein-ligand-docking/src/run_docking_pipeline.py \
  --config projects/04-protein-ligand-docking/configs/docking_config.json
```

### Run via Snakemake
```bash
snakemake --cores 2 -s projects/04-protein-ligand-docking/workflows/Snakefile
```

### Run Automated Pytest Suite
```bash
pytest projects/04-protein-ligand-docking/tests/ -v
```

---

## 📊 Generated Artifacts & Deliverables

- `data/target_ar_lbd.pdbqt`: Parameterized receptor model.
- `results/poses/*.pdbqt`: Multi-model PDBQTs containing all 9 conformational poses per compound.
- `results/docking_affinities.csv`: Ranked free binding energies ($\Delta G$).
- `results/interaction_contacts.csv`: Residue-level hydrogen bond and hydrophobic contact matrix.
- `results/docking_summary.json`: Parameter inventory, top affinities, and verification timestamps.
- `results/docking_analysis_report.html`: Interactive HTML summary dashboard with disclaimers.
- `figures/binding_affinity_comparison.png`: Predicted binding affinities comparison.
- `figures/interaction_heatmap.png`: Contact distances across active-site residues.
- `figures/pose_energy_spectrum.png`: Conformational mode rank vs. $\Delta G$ spectrum.
- `workflows/Dockerfile`: Isolated container definition for production deployments.
- `notebooks/04_molecular_docking_walkthrough.ipynb`: Interactive Jupyter walkthrough notebook.
