# Project 07: ESM-2 Deep Mutational Scanning & Zero-Shot Variant Fitness Prediction

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
[![PyTorch](https://img.shields.io/badge/deep%20learning-PyTorch%20%7C%20HuggingFace-EE4C2C.svg)](https://pytorch.org/)
[![Model: ESM-2](https://img.shields.io/badge/model-Meta%20ESM--2-blueviolet.svg)](https://github.com/facebookresearch/esm)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

Experimental deep mutational scanning (DMS) evaluates the functional impact of thousands of single amino-acid variants across a protein, but is labor-intensive, expensive, and limited by wet-lab assay scalability. Modern Protein Language Models (pLMs) trained on millions of evolutionary sequences capture deep structural, physical, and functional constraints directly in their attention representations.

This project implements a high-performance, GPU-accelerated **Variant Effect Prediction Engine** powered by Meta's **ESM-2** (Evolutionary Scale Modeling) transformer architecture. The pipeline performs zero-shot mutational fitness landscape scoring using masked-marginal log-likelihood ratios, benchmarks predictions against experimental DMS data from **MaveDB**, and visualizes mutational tolerance on 3D protein structures.

---

## Key Features

1. **Zero-Shot Masked Marginal Scoring (`src/esm_scorer.py`):**
   - Implements the masked marginal log-likelihood ratio (LLR) formulation:
     $$\Delta \log p(x_i) = \log p(x_i = y_{\text{mut}} \mid x_{\backslash i}) - \log p(x_i = y_{\text{wt}} \mid x_{\backslash i})$$
   - Computes fitness scores for all $19 \times L$ single amino-acid substitutions across sequence length $L$.
   - Supports batch tokenization, mixed-precision (`fp16`/`bf16`) inference, and flexible ESM-2 checkpoints (`esm2_t33_650M_UR50D`, `esm2_t12_35M_UR50D`, `esm2_t36_3B_UR50D`).

2. **DMS Experimental Benchmarking (`src/dms_evaluator.py`):**
   - Ingestion and normalization of gold-standard Deep Mutational Scanning datasets from **MaveDB** (e.g., *SARS-CoV-2 Spike RBD / ACE2 binding*, human *TP53* transactivation, and *EGFR* kinase inhibitor resistance).
   - Quantitative evaluation calculating Spearman rank correlation ($\rho$), Pearson correlation ($r$), and ROC-AUC / PR-AUC for classifying pathogenic vs. benign substitutions.
   - Direct comparison against legacy clinical pathogenicity scores (SIFT, PolyPhen-2, and CADD).

3. **Structural Mapping & Residue Tolerance (`src/structure_mapper.py`):**
   - Maps position-wise mutational tolerance (average mutational entropy / mean LLR) onto 3D experimental PDB structures or AlphaFold coordinates.
   - Computes residue solvent accessibility (SASA) and secondary structure to correlate pLM confidence with buried core vs. flexible surface loops.
   - Exports B-factor-replaced PDB files ready for ray-traced rendering in PyMOL or ChimeraX.

4. **Interactive Streamlit Web Dashboard (`src/app.py`):**
   - Clean, reactive web interface where users can upload any FASTA sequence or PDB ID.
   - Interactive 2D heatmap: residues ($1 \dots L$) vs. 20 amino acids with hoverable fitness values and classification flags.
   - Synchronized 3D structure viewer (powered by `py3Dmol`) dynamically color-coded by predicted mutational sensitivity.

---

## Directory Structure

```
07-protein-language-models/
├── configs/
│   └── esm_config.yaml           # Model checkpoint, batch size, device, and target genes
├── data/
│   ├── targets/                  # Target FASTA sequences & AlphaFold/PDB structures
│   └── mavedb_benchmarks/        # Ground-truth experimental DMS score matrices (CSV/TSV)
├── figures/                      # Heatmaps, correlation scatter plots & ROC curves
├── notebooks/
│   └── 07_esm2_variant_walkthrough.ipynb # Interactive tutorial & model inspection
├── results/
│   ├── dms_benchmark_metrics.json# Spearman rho, Pearson r, and ROC-AUC scores
│   ├── full_mutational_landscape.csv # Matrix of all 19 x L predicted LLRs
│   └── structural_tolerance_pdb/ # PDBs with B-factor replaced by mutational sensitivity
├── src/
│   ├── __init__.py
│   ├── esm_scorer.py             # PyTorch ESM-2 masked language inference engine
│   ├── dms_evaluator.py          # MaveDB parsing & statistical correlation metrics
│   ├── structure_mapper.py       # PDB B-factor injection & 3D coordinate mapping
│   ├── visualizer.py             # Publication heatmaps & benchmark plots
│   └── app.py                    # Interactive Streamlit dashboard
├── tests/
│   ├── test_esm_inference.py     # Unit tests for tokenization & score parity
│   ├── test_mavedb_parser.py     # Parsing tests for experimental score matrices
│   └── test_structure_mapper.py  # Validation of PDB B-factor substitution
└── workflows/
    └── Snakefile                 # Automated scoring and benchmark workflow
```

---

## Target Case Studies

- **Oncology Driver:** Human *TP53* DNA-binding domain (MaveDB: `urn:mavedb:00000068`).
- **Targeted Therapy Resistance:** Human *EGFR* kinase domain under Gefitinib/Osimertinib selection.
- **Viral Evolution:** SARS-CoV-2 Spike Receptor-Binding Domain (RBD) binding to human ACE2 (Starr et al., MaveDB: `urn:mavedb:00000030`).

---

## Quick Start

```bash
# 1. Run zero-shot mutational scoring on target sequence
python3 projects/07-protein-language-models/src/esm_scorer.py --fasta projects/07-protein-language-models/data/targets/tp53.fasta --output projects/07-protein-language-models/results/

# 2. Benchmark against experimental MaveDB data
python3 projects/07-protein-language-models/src/dms_evaluator.py --predicted projects/07-protein-language-models/results/full_mutational_landscape.csv --truth projects/07-protein-language-models/data/mavedb_benchmarks/tp53_dms.csv

# 3. Launch the interactive Streamlit dashboard
streamlit run projects/07-protein-language-models/src/app.py
```
