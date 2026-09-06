# Project 05: Transcriptomic Profiling & Differential Expression Analysis

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
[![Pipeline: Snakemake](https://img.shields.io/badge/workflow-Snakemake-brightgreen.svg)](https://snakemake.readthedocs.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This project implements an end-to-end transcriptomic analysis pipeline modeled on the **DESeq2** median-of-ratios normalization framework, negative binomial / generalized linear model differential expression testing, false discovery rate (FDR) control via Benjamini-Hochberg adjustment, and over-representation pathway enrichment analysis (ORA) using hypergeometric testing.

The dataset models transcriptional response in human muscle/androgen target cells (derived from **GEO accession GSE153664**) under control vs. dihydrotestosterone/androgen exposure across biological triplicates (3 Control vs. 3 Treated).

---

## Key Features

1. **Median-of-Ratios Normalization (`src/normalization.py`):**
   - Implements the DESeq2 size factor estimation algorithm: computes geometric pseudo-reference sample across all non-zero counts, calculates per-sample gene-to-reference ratios, and derives the median ratio as the library size factor.
   - Outputs both normalized count matrices and $\log_2(\text{norm\_count} + 1)$ transformed values for downstream variance stabilization.

2. **Differential Expression Testing (`src/differential_expression.py`):**
   - Estimates sample condition group means, pooled dispersion, and log2 fold change ($\log_2\text{FC}$).
   - Computes Wald $Z$-statistics and asymptotic two-sided $p$-values.
   - Adjusts for multiple testing using the **Benjamini-Hochberg (BH)** procedure to control the False Discovery Rate (FDR).
   - Flags statistically significant upregulated and downregulated transcripts according to user-defined thresholds (e.g., $|\log_2\text{FC}| \ge 1.0$, $\text{FDR} < 0.05$).

3. **Pathway Over-Representation Analysis (`src/pathway_enrichment.py`):**
   - Performs Fisher's exact / hypergeometric tests against defined gene sets (e.g., *Androgen Receptor Signaling Pathway*, *Skeletal Muscle Hypertrophy*, *Ribosome Biogenesis*, *Oxidative Phosphorylation*, *Fatty Acid Oxidation*).
   - Computes gene ratios, background ratios, hypergeometric $p$-values, and FDR values.

4. **Publication-Grade Visualizations (`src/visualizer.py`):**
   - **PCA Scatter Plot:** Principal component analysis on top 500 variable genes displaying sample clusters and variance explained.
   - **Volcano Plot:** Statistical significance ($-\log_{10}p$) vs. magnitude of change ($\log_2\text{FC}$) highlighting key induced/repressed biomarker genes (*KLK3*, *FKBP5*, *MYOD1*, *IGF1*, *SLC2A4*).
   - **Clustered Heatmap:** Clustered z-score expression profiles across samples for top differentially expressed genes.
   - **Pathway Enrichment Dotplot:** Gene ratio vs. $-\log_{10}(\text{FDR})$ with dot size mapped to gene count.

5. **Automated HTML Reporting (`shared/reporting/html_generator.py`):**
   - Generates an executive scientific report summarizing library sizes, DE transcript counts, pathway tables, and high-resolution figures.

---

## Directory Structure

```
05-rna-seq/
├── configs/
│   └── rnaseq_config.json        # Pipeline configuration & statistical thresholds
├── data/
│   ├── raw_counts_matrix.tsv     # GSE153664 raw transcript count matrix (547 genes x 6 samples)
│   └── sample_metadata.csv       # Sample metadata (sample_id, condition, replicate, sequencing depth)
├── figures/
│   ├── pca_plot.png              # PCA visualization
│   ├── volcano_plot.png          # Volcano plot of DE genes
│   ├── expression_heatmap.png    # Clustered expression heatmap
│   └── pathway_enrichment_dotplot.png # ORA dotplot
├── notebooks/
│   └── 05_rnaseq_de_walkthrough.ipynb # Interactive companion walkthrough
├── results/
│   ├── normalized_counts.csv     # DESeq2 normalized counts
│   ├── differential_expression_results.csv # Full statistical table
│   ├── significant_genes.csv     # Filtered significant hits
│   ├── pathway_enrichment_results.csv      # ORA pathway statistics
│   ├── rnaseq_summary.json       # Machine-readable JSON summary
│   └── rnaseq_analysis_report.html # Comprehensive HTML report
├── src/
│   ├── __init__.py
│   ├── normalization.py          # Median-of-ratios algorithm
│   ├── differential_expression.py # Wald test & BH FDR adjustment
│   ├── pathway_enrichment.py     # Hypergeometric ORA
│   ├── visualizer.py             # Publication graphics
│   └── run_rnaseq_pipeline.py    # Master CLI pipeline runner
├── tests/
│   ├── test_normalization.py     # Unit tests for size factors & scaling
│   ├── test_de.py                # Unit tests for Wald test & BH procedure
│   └── test_enrichment.py        # Unit tests for hypergeometric test
└── workflows/
    └── Snakefile                 # Reproducible Snakemake workflow
```

---

## Quick Start

### 1. Run the Pipeline

```bash
python3 projects/05-rna-seq/src/run_rnaseq_pipeline.py --config projects/05-rna-seq/configs/rnaseq_config.json
```

### 2. Run Tests

```bash
pytest projects/05-rna-seq/tests/ -v
```

### 3. Open the Interactive Notebook

```bash
jupyter lab projects/05-rna-seq/notebooks/05_rnaseq_de_walkthrough.ipynb
```

---

## Scientific Context & Disclaimer

> [!NOTE]
> This analysis is conducted exclusively for computational biology methodology demonstration and scientific transcriptomic exploration. All gene sets and expression matrices reflect public research benchmarks (GEO GSE153664) and do not constitute clinical guidance or endocrine therapy protocols.
