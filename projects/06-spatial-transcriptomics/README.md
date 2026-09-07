# Project 06: Spatial Transcriptomics & Single-Cell Deconvolution

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
[![Scanpy & Squidpy](https://img.shields.io/badge/single--cell-Scanpy%20%7C%20Squidpy-brightgreen.svg)](https://scanpy.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

Single-cell RNA sequencing (scRNA-seq) offers single-cell resolution but destroys spatial tissue context during tissue dissociation. Conversely, commercial spatial transcriptomics technologies (such as 10x Genomics Visium) preserve spatial architecture and histology but capture 1–10 cells per 55µm spot.

This project delivers an end-to-end, production-grade **Spatial & Single-Cell Deconvolution Pipeline** to dissect the tumor microenvironment (TME) in human cancer tissue. By mapping high-resolution scRNA-seq reference signatures onto spatially resolved transcriptomic slides, this pipeline infers continuous cell-type abundances, constructs spatial tissue graphs, and models paracrine ligand-receptor communications across physical tumor-stroma boundaries.

---

## Key Features

1. **Single-Cell Reference Quality Control & Annotation (`src/qc_reference.py`):**
   - Automated QC filtering on counts per cell, detected genes, mitochondrial percentage, and doublet identification via `Scrublet`.
   - Normalization, highly variable gene (HVG) selection, PCA, and batch integration using `Harmony` or deep generative modeling (`scVI`).
   - Leiden/Louvain clustering and automated cluster annotation via canonical cell-type markers (e.g., CD8+ T cells, regulatory T cells, cancer-associated fibroblasts, malignant epithelial cells).

2. **Spatial Data Preprocessing & Histology Alignment (`src/spatial_qc.py`):**
   - High-resolution H&E image processing and spot coordinate alignment using `Squidpy` and `AnnData`.
   - Spot-level filtering (library size, spot gene counts, artifact masking).
   - Spatially variable gene (SVG) discovery using spatial autocorrelation metrics (Moran's $I$ and Geary's $C$).

3. **Probabilistic Cell-Type Deconvolution (`src/deconvolution.py`):**
   - Implements Bayesian/probabilistic reference-based deconvolution (integrating principles from `Cell2location` and negative binomial regression).
   - Maps cell-type abundance distributions onto every Visium capture spot, estimating absolute cell densities rather than relative proportions alone.

4. **Spatial Cell-Cell Communication (`src/spatial_interactions.py`):**
   - Spatial neighborhood graph construction (Delaunay triangulation / radial distance graphs).
   - Co-localization analysis evaluating immune-cell exclusion and stromal trapping.
   - Paracrine signaling inference linking spatially adjacent ligand-receptor pairs (e.g., *PD-L1 / PD-1*, *CXCL12 / CXCR4*, *TGFB1 / TGFBR2*) within defined physical micro-domains.

5. **Interactive Visualization & Spatial Dashboards (`src/visualizer.py`):**
   - Histology overlays rendering deconvolved cell densities directly onto high-resolution H&E tissue images.
   - Spatial feature expression maps, spatial Moran's $I$ rank plots, and interactive cell-cell communication chord diagrams.

---

## Directory Structure

```
06-spatial-transcriptomics/
├── configs/
│   └── spatial_config.yaml       # Quality thresholds, model hyperparameters & gene sets
├── data/
│   ├── reference_sc/             # Matched scRNA-seq h5ad reference count matrix
│   └── spatial_visium/           # 10x Visium spatial folder (filtered_feature_bc_matrix.h5, spatial/)
├── figures/                      # Deconvolution histology overlays & interaction plots
├── notebooks/
│   └── 06_spatial_deconvolution_walkthrough.ipynb # Interactive end-to-end tutorial
├── results/
│   ├── cell_abundances_per_spot.csv  # Deconvolved cell density matrix
│   ├── spatially_variable_genes.csv  # Moran's I and p-values
│   └── spatial_interactions.csv      # Prioritized ligand-receptor spatial interactions
├── src/
│   ├── __init__.py
│   ├── qc_reference.py           # scRNA-seq filtering, clustering, and signature extraction
│   ├── spatial_qc.py             # Visium spot filtering & spatial autocorrelation
│   ├── deconvolution.py          # Probabilistic spot deconvolution engine
│   ├── spatial_interactions.py   # Neighborhood graph & ligand-receptor analysis
│   └── visualizer.py             # Spatial plotting & H&E overlay generation
├── tests/
│   ├── test_qc.py                # Unit tests for filtering thresholds and metrics
│   ├── test_deconvolution.py     # Verification of deconvolution output dimensions and non-negativity
│   └── test_spatial_graph.py     # Verification of spatial coordinate adjacency graphs
└── workflows/
    └── Snakefile                 # Reproducible workflow pipeline
```

---

## Target Datasets

- **Spatial Platform:** 10x Genomics Visium (Human Breast Cancer or Colorectal Carcinoma).
- **Single-Cell Reference:** Matched human tumor atlas (e.g., GSE176078 or GEO/Zenodo public Visium benchmarks).

---

## Quick Start

```bash
# 1. Run the deconvolution workflow
python3 projects/06-spatial-transcriptomics/src/run_spatial_pipeline.py --config projects/06-spatial-transcriptomics/configs/spatial_config.yaml

# 2. Run test suite
pytest projects/06-spatial-transcriptomics/tests/ -v

# 3. Explore interactive walkthrough
jupyter lab projects/06-spatial-transcriptomics/notebooks/06_spatial_deconvolution_walkthrough.ipynb
```
