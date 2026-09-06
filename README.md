# 🧬 Computational Biology & Bioinformatics Engineering Portfolio

### *Software Engineering × Computational Biology × Bioinformatics*

[![CI](https://github.com/example/bioinformatics-portfolio/actions/workflows/ci.yml/badge.svg)](https://github.com/example/bioinformatics-portfolio/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![Code Style: Ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **"I combine software engineering with computational biology to build reproducible tools and workflows for biological-data analysis."**

---

## 🎯 Engineering & Biological Positioning

This repository showcases production-grade scientific software engineering applied to fundamental challenges in molecular biology, structural pharmacology, evolutionary genomics, and transcriptomics. Rather than disconnected tutorials or shallow notebooks, this portfolio demonstrates **robust system architecture**, **deterministic algorithms**, **clean test suites**, and **reproducible workflows** utilizing real public biological data from NCBI, UniProt, RCSB PDB, AlphaFold DB, and GEO.

```
                     BIOLOGICAL & COMPUTATIONAL PROGRESSION
  
  01. SEQUENCE        02. EVOLUTION       03. STRUCTURE       04. INTERACTION     05. TRANSCRIPTOMICS
  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐    ┌──────────────┐    ┌─────────────────┐
  │ FASTA Parser │    │ Multi-Align  │    │ AlphaFold v2 │    │ AutoDock     │    │ RNA-Seq Counts  │
  │ Physicochem  ├───►│ Orthologs    ├───►│ pLDDT / PAE  ├───►│ Vina Docking ├───►│ Normalization   │
  │ BLAST+ CLI   │    │ Phylogeny    │    │ Pocket Cleft │    │ Interaction  │    │ DE & Pathways   │
  └──────────────┘    └──────────────┘    └──────────────┘    └──────────────┘    └─────────────────┘
```

---

## 📂 Flagship Projects Matrix

| Project | Core Domain | Primary Tools & Algorithms | Key Engineering Deliverables | Real Dataset / Accession |
| :--- | :--- | :--- | :--- | :--- |
| **[01. Sequence Homology](projects/01-sequence-homology)** | Sequence Analytics & Homology | Python 3.11, Biopython, Typer, BLAST+, Rich | Production CLI (`bioseq`), physicochemical profiler, batch BLAST engine, HTML report | Human AR (`UniProt: P10275`), Insulin (`P01308`) |
| **[02. AR Phylogenetics](projects/02-androgen-receptor-phylogenetics)** | Evolutionary Genomics | MAFFT, MUSCLE, IQ-TREE, Biopython Phylo | Ortholog alignment pipeline, Shannon entropy scoring, Newick tree visualizer | Vertebrate AR Panel (`NCBI: NM_000044`, 12 Species) |
| **[03. AlphaFold Structure](projects/03-alphafold-structure-analysis)** | Structural Bioinformatics | AlphaFold v2, PDB, DSSP, PyMOL, ChimeraX | pLDDT/PAE residue profiler, secondary structure calculator, ray-trace render scripts | Human AR LBD (`PDB: 1E3G`, `AF-P10275-F1`) |
| **[04. Molecular Docking](projects/04-protein-ligand-docking)** | Computational Pharmacology | AutoDock Vina, RDKit, PLIP, Snakemake | Parameterized docking pipeline, MMFF94 ligand prep, interaction fingerprinting | AR-LBD + DHT / Enzalutamide (`PubChem: 10635`) |
| **[05. RNA-Seq Pipeline](projects/05-rna-seq)** | Functional Transcriptomics | DESeq2 / Negative Binomial GLM, ORA, Snakemake | FastQC QC synthesis, Median-of-Ratios norm, Wald test DE, Volcano/PCA, GO/KEGG | Hormone Response Muscle Study (`GEO: GSE153664`) |

---

## 🏛️ System Architecture

```
bioinformatics-portfolio/
├── shared/                         # Reusable core engineering libraries
│   ├── io/                         # FASTA, PDB, PDBQT, Tabular parsers
│   ├── logging/                    # Structured logging with Rich formatters
│   ├── visualization/              # Publication themes (300 DPI, vector export)
│   ├── reporting/                  # Standalone Jinja2 HTML report generator
│   └── utils/                      # Rate-limiting, checksums, physical constants
│
├── projects/                       # 5 Flagship computational biology projects
│   ├── 01-sequence-homology/       # bioseq CLI tool and BLAST orchestration
│   ├── 02-androgen-receptor-phylogenetics/ # Comparative MSA and tree inference
│   ├── 03-alphafold-structure-analysis/    # 3D structure and pLDDT/PAE metrics
│   ├── 04-protein-ligand-docking/  # AutoDock Vina & PLIP interaction workflow
│   └── 05-rna-seq/                 # End-to-end RNA-seq DE & pathway analysis
│
├── docs/                           # Architecture, reproducibility & services
│   ├── reproducibility.md          # Provenance tracking & deterministic guidelines
│   ├── tools.md                    # Algorithmic & mathematical inventory
│   ├── architecture.md             # Clean architecture specifications
│   └── freelancing-services.md     # Client service mapping & ethical boundaries
│
├── tests/                          # Root integration & shared tests
├── Makefile                        # Developer & workflow orchestration
└── pyproject.toml                  # Pinned dependencies & build configuration
```

---

## 🛠️ Technology Stack

- **Core & Scientific Computing**: Python 3.11+, NumPy, SciPy, Pandas, Statsmodels.
- **Bioinformatics Libraries**: Biopython, RDKit, PyMOL automation scripts.
- **Workflow & CLI Engineering**: Typer, Rich, Snakemake, Jinja2.
- **Visualization**: Matplotlib, Seaborn, SVG/HTML dynamic reporting.
- **Quality & CI/CD**: Pytest, Pytest-Cov, Ruff, Mypy, GitHub Actions.

---

## ⚡ Quick Start

### 1. Installation

```bash
# Clone the repository
git clone https://github.com/aashiq-parinda/bioinformatics-portfolio.git
cd bioinformatics-portfolio

# Create and activate a virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install all dependencies and CLI in editable mode
make install-dev
```

### 2. Run Test Suite & Linters

```bash
# Run all tests with pytest
make test

# Run Ruff linting and formatting verification
make check
```

### 3. Run Flagship Pipelines

```bash
# Project 01: Analyze protein FASTA and run BLAST pipeline
bioseq analyze projects/01-sequence-homology/examples/human_androgen_receptor.fasta --output projects/01-sequence-homology/results/

# Project 02: Run Androgen Receptor Phylogenetic Pipeline
make run-p02

# Project 03: Run AlphaFold Structure & Confidence Assessment
make run-p03

# Project 04: Execute Molecular Docking Pipeline
make run-p04

# Project 05: Execute End-to-End RNA-Seq Analysis
make run-p05

# Or run all flagship pipelines sequentially
make run-all
```

---

## 🔬 Scientific Safety & Ethical Disclaimer

1. **Computational Research Only**: The models, docking simulations, and phylogenetic analyses in this repository are intended solely for computational biology research, tool development, and educational evaluation.
2. **No Clinical or Medical Advice**: None of the code or documentation provides medical diagnosis, therapy recommendations, or anabolic cycle/dosing instructions.
3. **Hypothesis Generation**: Computational molecular docking affinities are theoretical hypotheses requiring experimental wet-lab validation (e.g., surface plasmon resonance or isothermal titration calorimetry).

---

## 💼 Freelance & Collaboration Services

Looking for robust computational biology tools, automated bioinformatics pipelines, or custom data analysis workflows? See [docs/freelancing-services.md](docs/freelancing-services.md) for available client engagements.

## 📄 License & Citation

This project is licensed under the MIT License - see [LICENSE](LICENSE) for details. To cite this repository in academic or computational work, see [CITATION.cff](CITATION.cff).
