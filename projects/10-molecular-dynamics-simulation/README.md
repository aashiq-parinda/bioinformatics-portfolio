# Project 10: Explicit-Solvent Molecular Dynamics & Binding Free Energy (OpenMM & MM-PBSA)

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
[![OpenMM](https://img.shields.io/badge/simulation-OpenMM%20%7C%20CUDA-green.svg)](https://openmm.org/)
[![MDAnalysis](https://img.shields.io/badge/analysis-MDAnalysis-red.svg)](https://www.mdanalysis.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

Static molecular docking algorithms (such as AutoDock Vina in Project 04) provide fast conformational screening but suffer from high false-positive rates due to rigid-receptor approximations and empirical scoring functions. In modern structure-based drug design (SBDD), putative hits must be validated using all-atom, explicit-solvent **Molecular Dynamics (MD)** simulations to capture receptor plasticity, examine solvent-mediated hydrogen-bond networks, assess ligand residence stability, and compute rigorous thermodynamic binding free energies.

This project implements an automated, production-grade **Molecular Dynamics & Free Energy Assessment Pipeline** powered by **OpenMM** and **MDAnalysis**. The pipeline automates force field parameterization of novel small-molecule ligands, constructs solvated periodic boundary boxes, executes multi-stage equilibration and production runs, extracts biophysical trajectory metrics, and computes binding free energy ($\Delta G_{\text{bind}}$) using the **MM-PBSA / MM-GBSA** methodology.

---

## Biophysical Protocol & Simulation Workflow

```
       DOCKED COMPLEX (PDB + SDF)
                   │
       ┌───────────┴───────────┐
       ▼                       ▼
  PROTEIN SETUP           LIGAND PARAMETERIZATION
  (AMBER ff14SB)          (OpenFF 2.1.0 / GAFF2 + AM1-BCC Charges)
       │                       │
       └───────────┬───────────┘
                   ▼
       SOLVATION & IONIZATION
       (TIP3P water box, 10Å padding, 0.15 M NaCl)
                   │
                   ▼
       ENERGY MINIMIZATION
       (Conjugate gradient / L-BFGS to tolerance < 10 kJ/mol/nm)
                   │
                   ▼
       NVT THERMAL EQUILIBRATION
       (Heating to 300 K with harmonic position restraints)
                   │
                   ▼
       NPT PRESSURE EQUILIBRATION
       (1.0 atm Monte Carlo barostat with relaxed restraints)
                   │
                   ▼
       PRODUCTION DYNAMICS (OpenMM GPU Engine)
       (20 - 100 ns unrestrained NPT, 2 fs timestep, Langevin dynamics)
                   │
       ┌───────────┴────────────────────────┐
       ▼                                    ▼
  TRAJECTORY ANALYSIS               BINDING FREE ENERGY
  - Protein/Ligand RMSD             (MM-PBSA / MM-GBSA)
  - Residue RMSF                    - Electrostatic & vdW terms
  - H-Bond Persistence Fraction     - Polar & Non-polar solvation
  - Radius of Gyration (Rg)         - Net ΔG_bind calculation
```

---

## Key Features

1. **Automated Small Molecule & Protein Parameterization (`src/system_builder.py`):**
   - Integrates the **Open Force Field Toolkit (`OpenFF`)** and `RDKit` to assign AM1-BCC partial charges and parameters to small molecules without manual intervention.
   - Cleans receptor structures, models missing loops/atoms via `PDBFixer`, and applies `AMBER14` protein force fields.
   - Solvates systems in rhombic dodecahedron or cubic periodic boundary condition (PBC) boxes with explicit TIP3P water and neutralizes with physiological ionic strength (0.15 M NaCl).

2. **GPU-Accelerated Simulation Engine (`src/md_engine.py`):**
   - Fully scriptable OpenMM engine utilizing OpenCL/CUDA acceleration with Particle Mesh Ewald (PME) electrostatics and 9Å non-bonded cutoffs.
   - Automated four-stage protocol: Energy Minimization $\to$ NVT Heating $\to$ NPT Density Equilibration $\to$ Unrestrained Production.
   - Regular trajectory checkpointing (`.dcd` / `.xtc`) and state logging (potential energy, kinetic energy, temperature, density).

3. **Biophysical Trajectory Analytics (`src/trajectory_analyzer.py`):**
   - Streamed parsing of multi-nanosecond trajectories via `MDAnalysis`.
   - Calculates time-dependent Root-Mean-Square Deviation (RMSD) for protein backbone, binding pocket residues, and ligand heavy atoms.
   - Computes per-residue Root-Mean-Square Fluctuation (RMSF) to identify induced stabilization.
   - Evaluates hydrogen bond formation frequency, donor-acceptor distances, and persistence percentages across the simulation timeline.

4. **Thermodynamic Free Energy Calculation (`src/mmpbsa_calculator.py`):**
   - End-state free energy calculation:
     $$\Delta G_{\text{bind}} = \langle G_{\text{complex}} \rangle - \langle G_{\text{protein}} \rangle - \langle G_{\text{ligand}} \rangle$$
   - Solves the Poisson-Boltzmann (PB) or Generalized Born (GB) equations for polar solvation and surface area (SASA) for nonpolar solvation.
   - Provides per-residue energy decomposition identifying critical hotspot residues driving binding affinity.

5. **Visual Reporting & Trajectory Rendering (`src/visualizer.py`):**
   - Publication-quality plots: RMSD/RMSF curves, 2D RMSD cross-correlation heatmaps, and H-bond occupancy bar charts.
   - Trajectory alignment and automated PyMOL rendering scripts to produce smooth rotation animations and contact maps.

---

## Directory Structure

```
10-molecular-dynamics-simulation/
├── configs/
│   └── simulation_config.yaml    # Temperature, pressure, timestep, duration & force field specs
├── data/
│   ├── complexes/                # Input receptor-ligand complexes (PDB + SDF)
│   └── topology/                 # Generated solvated topology and parameter files
├── figures/                      # RMSD, RMSF, H-bond, and free energy plots
├── notebooks/
│   └── 10_molecular_dynamics_walkthrough.ipynb # Interactive simulation inspection
├── results/
│   ├── trajectories/             # DCD/XTC coordinate trajectories & logs
│   ├── analysis_metrics.csv      # Time-series RMSD, RMSF, and contact data
│   └── mmpbsa_free_energy.json   # Decomposed free energy terms (ΔE_MM, ΔG_solv, ΔG_bind)
├── src/
│   ├── __init__.py
│   ├── system_builder.py         # OpenFF ligand prep, PDBFixer, and box solvation
│   ├── md_engine.py              # OpenMM minimization, equilibration, and production
│   ├── trajectory_analyzer.py    # MDAnalysis RMSD, RMSF, and H-bond profiling
│   ├── mmpbsa_calculator.py      # MM-PBSA / MM-GBSA free energy calculator
│   └── visualizer.py             # Trajectory plotting and animation generator
├── tests/
│   ├── test_parameterization.py  # Unit tests for small-molecule force field assignment
│   ├── test_system_builder.py    # Unit tests for box dimension and charge neutrality
│   └── test_trajectory_analysis.py # Unit tests for RMSD and RMSF computation
└── workflows/
    └── Snakefile                 # End-to-end automated MD pipeline
```

---

## Target Case Study

- **Receptor:** Human Androgen Receptor Ligand-Binding Domain (AR-LBD, `PDB: 1E3G`).
- **Ligands:** Endogenous agonist Dihydrotestosterone (DHT) vs. second-generation antiandrogen Enzalutamide vs. resistance-inducing mutants (e.g., F877L).

---

## Quick Start

```bash
# 1. Parameterize and solvate target protein-ligand complex
python3 projects/10-molecular-dynamics-simulation/src/system_builder.py \
    --protein projects/10-molecular-dynamics-simulation/data/complexes/ar_lbd.pdb \
    --ligand projects/10-molecular-dynamics-simulation/data/complexes/enzalutamide.sdf \
    --output projects/10-molecular-dynamics-simulation/data/topology/

# 2. Run OpenMM simulation pipeline
python3 projects/10-molecular-dynamics-simulation/src/md_engine.py \
    --config projects/10-molecular-dynamics-simulation/configs/simulation_config.yaml

# 3. Analyze trajectory and compute MM-PBSA binding free energy
python3 projects/10-molecular-dynamics-simulation/src/trajectory_analyzer.py \
    --topology projects/10-molecular-dynamics-simulation/data/topology/system.pdb \
    --trajectory projects/10-molecular-dynamics-simulation/results/trajectories/production.dcd
```
