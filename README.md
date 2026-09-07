# 🧬 Computational Biology & Bioinformatics Engineering Portfolio

### *Software Engineering × Computational Biology × Bioinformatics*

[![CI](https://github.com/example/bioinformatics-portfolio/actions/workflows/ci.yml/badge.svg)](https://github.com/example/bioinformatics-portfolio/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![Code Style: Ruff](https://img.shields.io/badge/code%20style-ruff-000000.svg)](https://github.com/astral-sh/ruff)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **"I combine software engineering with computational biology to build reproducible tools and workflows for biological-data analysis."**

---

## 🎯 Engineering & Biological Positioning

This repository showcases production-grade scientific software engineering applied to fundamental challenges in molecular biology, structural pharmacology, evolutionary genomics, transcriptomics, generative AI, and clinical bioinformatics. Rather than disconnected tutorials or shallow notebooks, this portfolio demonstrates **robust system architecture**, **deterministic algorithms**, **clean test suites**, and **reproducible production workflows** utilizing real public biological data from NCBI, UniProt, RCSB PDB, AlphaFold DB, GEO, and MaveDB.

```
                     BIOLOGICAL & COMPUTATIONAL PROGRESSION

  [PART I: FOUNDATIONAL GENOMICS & STRUCTURAL PHARMACOLOGY]
  01. SEQUENCE        02. EVOLUTION       03. STRUCTURE       04. DOCKING         05. RNA-SEQ
  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐    ┌──────────────┐    ┌─────────────────┐
  │ FASTA Parser │    │ Multi-Align  │    │ AlphaFold v2 │    │ AutoDock     │    │ RNA-Seq Counts  │
  │ Physicochem  ├───►│ Orthologs    ├───►│ pLDDT / PAE  ├───►│ Vina Docking ├───►│ Normalization   │
  │ BLAST+ CLI   │    │ Phylogeny    │    │ Pocket Cleft │    │ Interaction  │    │ DE & Pathways   │
  └──────────────┘    └──────────────┘    └──────────────┘    └──────────────┘    └─────────────────┘
                                                                                            │
  [PART II: MODERN MODALITIES, GENERATIVE AI & PRODUCTION ENGINEERING]                      ▼
  10. DYNAMICS        09. CLINICAL        08. SURVEILLANCE    07. AI PROTEIN      06. SPATIAL OMICS
  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐    ┌──────────────┐    ┌─────────────────┐
  │ OpenMM GPU   │    │ AMP / ASCO   │    │ Nextflow DSL2│    │ Meta ESM-2   │    │ 10x Visium + sc │
  │ MM-PBSA Free │◄───┤ VCF Somatic  │◄───┤ AMR Profiling│◄───┤ Zero-Shot DMS│◄───┤ Cell2location   │
  │ Energy & RMSD│    │ CIViC Action │    │ MultiQC Cloud│    │ MaveDB Bench │    │ Ligand-Receptor │
  └──────────────┘    └──────────────┘    └──────────────┘    └──────────────┘    └─────────────────┘
```

---

## 📂 Flagship Projects Matrix

### Part I: Foundational Genomics & Structural Pharmacology

| Project | Core Domain | Primary Tools & Algorithms | Key Engineering Deliverables | Real Dataset / Accession |
| :--- | :--- | :--- | :--- | :--- |
| **[01. Sequence Homology](projects/01-sequence-homology)** | Sequence Analytics & Homology | Python 3.11, Biopython, Typer, BLAST+, Rich | Production CLI (`bioseq`), physicochemical profiler, batch BLAST engine, HTML report | Human AR (`UniProt: P10275`), Insulin (`P01308`) |
| **[02. AR Phylogenetics](projects/02-androgen-receptor-phylogenetics)** | Evolutionary Genomics | MAFFT, MUSCLE, IQ-TREE, Biopython Phylo | Ortholog alignment pipeline, Shannon entropy scoring, Newick tree visualizer | Vertebrate AR Panel (`NCBI: NM_000044`, 12 Species) |
| **[03. AlphaFold Structure](projects/03-alphafold-structure-analysis)** | Structural Bioinformatics | AlphaFold v2, PDB, DSSP, PyMOL, ChimeraX | pLDDT/PAE residue profiler, secondary structure calculator, ray-trace render scripts | Human AR LBD (`PDB: 1E3G`, `AF-P10275-F1`) |
| **[04. Molecular Docking](projects/04-protein-ligand-docking)** | Computational Pharmacology | AutoDock Vina, RDKit, PLIP, Snakemake | Parameterized docking pipeline, MMFF94 ligand prep, interaction fingerprinting | AR-LBD + DHT / Enzalutamide (`PubChem: 10635`) |
| **[05. RNA-Seq Pipeline](projects/05-rna-seq)** | Functional Transcriptomics | DESeq2 / Negative Binomial GLM, ORA, Snakemake | FastQC QC synthesis, Median-of-Ratios norm, Wald test DE, Volcano/PCA, GO/KEGG | Hormone Response Muscle Study (`GEO: GSE153664`) |

### Part II: Modern Modalities, Generative AI & Production Systems

| Project | Core Domain | Primary Tools & Algorithms | Key Engineering Deliverables | Real Dataset / Accession |
| :--- | :--- | :--- | :--- | :--- |
| **[06. Spatial Transcriptomics](projects/06-spatial-transcriptomics)** | Single-Cell & Spatial Omics | Scanpy, Squidpy, AnnData, scVI, Cell2location | Tumor microenvironment deconvolution, Moran's $I$ spatial autocorrelation, paracrine signaling | 10x Visium Breast/CRC + matched scRNA-seq (`GSE176078`) |
| **[07. Protein Language Models](projects/07-protein-language-models)** | Generative AI & Protein Engineering | PyTorch, Meta ESM-2, HuggingFace, Streamlit | Zero-shot masked marginal LLR fitness scoring, MaveDB DMS validation, interactive 3D mutational dashboard | Human *TP53* (`urn:mavedb:00000068`), SARS-CoV-2 RBD |
| **[08. Nextflow Surveillance](projects/08-nextflow-pathogen-surveillance)** | Bioinformatics DevOps & Cloud Pipelines | Nextflow DSL2, Docker/Singularity, MultiQC | nf-core compliant pathogen surveillance pipeline, AMR gene profiling, automated CI/CD & cloud scaling | Viral/Microbial WGS (Illumina & Oxford Nanopore) |
| **[09. Clinical Interpretation](projects/09-clinical-variant-interpretation)** | Clinical Genomics & Precision Oncology | Python, cyvcf2, Ensembl VEP, CIViC REST API | Automated AMP/ASCO/CAP 4-tier somatic curation, drug matching engine, molecular pathology PDF report | Somatic NGS Tumor-Normal Panels (`MSK-IMPACT` / TCGA) |
| **[10. Molecular Dynamics](projects/10-molecular-dynamics-simulation)** | Computational Biophysics & Drug Discovery | OpenMM, MDAnalysis, OpenFF, AMBER ff14SB | Explicit-solvent GPU MD engine, RMSD/RMSF & H-bond stability, MM-PBSA binding free energy ($\Delta G_{\text{bind}}$) | AR-LBD + Enzalutamide / Resistance Mutants (`PDB: 1E3G`) |

---

## 🏛️ System Architecture

```
bioinformatics-portfolio/
├── shared/                             # Reusable core engineering libraries
│   ├── io/                             # FASTA, PDB, PDBQT, Tabular parsers
│   ├── logging/                        # Structured logging with Rich formatters
│   ├── visualization/                  # Publication themes (300 DPI, vector export)
│   ├── reporting/                      # Standalone Jinja2 HTML report generator
│   └── utils/                          # Rate-limiting, checksums, physical constants
│
├── projects/                           # 10 Flagship computational biology projects
│   ├── 01-sequence-homology/           # bioseq CLI tool and BLAST orchestration
│   ├── 02-androgen-receptor-phylogenetics/ # Comparative MSA and tree inference
│   ├── 03-alphafold-structure-analysis/    # 3D structure and pLDDT/PAE metrics
│   ├── 04-protein-ligand-docking/      # AutoDock Vina & PLIP interaction workflow
│   ├── 05-rna-seq/                     # End-to-end RNA-seq DE & pathway analysis
│   ├── 06-spatial-transcriptomics/     # 10x Visium & scRNA-seq tumor deconvolution
│   ├── 07-protein-language-models/     # ESM-2 zero-shot DMS fitness prediction & app
│   ├── 08-nextflow-pathogen-surveillance/ # Nextflow DSL2 pathogen WGS & AMR pipeline
│   ├── 09-clinical-variant-interpretation/# AMP/ASCO/CAP oncology variant curation engine
│   └── 10-molecular-dynamics-simulation/  # OpenMM explicit-solvent MD & MM-PBSA
│
├── docs/                               # Architecture, reproducibility & services
├── tests/                              # Root integration & shared tests
├── Makefile                            # Developer & workflow orchestration
└── pyproject.toml                      # Pinned dependencies & build configuration
```

---

## 🛠️ Technology Stack

- **Deep Learning & Generative Biology**: PyTorch, Hugging Face Transformers, Meta ESM-2, scVI-tools.
- **Single-Cell & Spatial Omics**: Scanpy, Squidpy, AnnData, Cell2location.
- **Biophysics & Molecular Modeling**: OpenMM, MDAnalysis, OpenFF, RDKit, Biopython.
- **Pipeline Orchestration & DevOps**: Nextflow (DSL2), Snakemake, Docker, Singularity, MultiQC, GitHub Actions CI.
- **Clinical Genomics**: cyvcf2, pysam, Ensembl VEP, CIViC API, ReportLab.
- **Web & Interactive Analytics**: Streamlit, py3Dmol, Rich, Typer, Jinja2.

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
