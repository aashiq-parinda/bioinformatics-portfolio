# Project 03: Protein Structure & AlphaFold Analysis

[![Dataset: AlphaFold DB](https://img.shields.io/badge/AlphaFold%20DB-AF--P10275--F1-blue.svg)](data/af_p10275_human_ar.pdb)
[![Dataset: RCSB PDB](https://img.shields.io/badge/RCSB%20PDB-1E3G-darkgreen.svg)](data/1e3g_human_ar_lbd.pdb)
[![Tests: Pytest](https://img.shields.io/badge/tests-passing-brightgreen.svg)](tests/)

A computational structural biology investigation assessing AlphaFold v2 prediction metrics, per-residue confidence (pLDDT) distributions, secondary structure assignments, and experimental crystallographic validation against the 2.4 Å X-ray crystal structure of the Human Androgen Receptor (PDB: `1E3G`).

---

## 🧬 Scientific Rationale & Structural Biology Questions

Nuclear receptor structures demonstrate extreme modularity:
1. **Intrinsically Disordered N-Terminal Domain (NTD, 1-555)**: Does AlphaFold capture its conformational plasticity by reporting low pLDDT ($<50$)?
2. **Folded Functional Cores (DBD 556-623, LBD 676-919)**: Does AlphaFold achieve high-confidence prediction ($>90$) for the closed 12-$\alpha$-helix steroid binding sandwich?
3. **Active Pocket Fidelity**: Do the key steroid-coordinating residues (Asn705, Thr877, Phe764, Met745) in the experimental DHT-bound cavity align with AlphaFold's modeled side chains?

---

## 🔬 Benchmark Structural Datasets

- **`data/af_p10275_human_ar.pdb`**: Full-length 919-residue human androgen receptor from AlphaFold DB (`AF-P10275-F1-model_v6.pdb`).
- **`data/1e3g_human_ar_lbd.pdb`**: High-resolution X-ray crystal structure of Human AR LBD complexed with $5\alpha$-dihydrotestosterone (DHT) from RCSB PDB (`DOI: 10.2210/pdb1e3g/pdb`).

---

## 🛠️ Execution & Commands

### 1. Run Complete Analysis Pipeline
```bash
python3 projects/03-alphafold-structure-analysis/src/run_analysis.py \
  --config projects/03-alphafold-structure-analysis/configs/analysis_config.json
```

### 2. Run Pytest Suite
```bash
pytest projects/03-alphafold-structure-analysis/tests/ -v
```

### 3. PyMOL & ChimeraX Automation
```bash
# Launch in PyMOL
pymol projects/03-alphafold-structure-analysis/scripts/visualize_ar_structure.pml

# Launch in ChimeraX
chimerax projects/03-alphafold-structure-analysis/scripts/visualize_ar_structure.cxc
```

---

## 📊 Pipeline Deliverables & Generated Artifacts

- `results/per_residue_plddt.csv`: Per-residue pLDDT score and confidence tier classification.
- `results/secondary_structure_assignments.csv`: Dihedral angle ($\phi, \psi$) and secondary structure states.
- `results/pocket_residues.csv`: Residues within 4.5 Å of ligand DHT with AlphaFold pLDDT scores.
- `results/structure_summary.json`: Global metrics, RMSD vs 1E3G, and domain confidence values.
- `results/alphafold_structure_report.html`: Interactive HTML report with data tables and figures.
- `figures/plddt_residue_profile.png`: Per-residue pLDDT with domain demarcations.
- `figures/confidence_distribution.png`: Confidence category proportions bar chart.
- `figures/pocket_residues_plddt.png`: Confidence of catalytic binding pocket residues.
- `notebooks/03_alphafold_structure_walkthrough.ipynb`: Interactive Jupyter walkthrough notebook.
