# Software Architecture & Engineering Standards

## 1. Architectural Philosophy

This repository is architected following the clean architecture pattern adapted for scientific workflows and high-performance computational biology. The design cleanly separates **data ingestion**, **analytical core logic**, **scientific visualization**, and **user-facing CLI/orchestration interfaces**.

```
┌────────────────────────────────────────────────────────────────────────┐
│                        PRESENTATION & CLI                             │
│       Typer CLI (bioseq)  •  Snakemake Workflows  •  HTML Reports      │
└───────────────────────────────────┬────────────────────────────────────┘
                                    │
┌───────────────────────────────────▼────────────────────────────────────┐
│                       DOMAIN PIPELINES & WORKFLOWS                     │
│  01-Homology   02-Phylogenetics   03-AlphaFold   04-Docking   05-RNA   │
└───────────────────────────────────┬────────────────────────────────────┘
                                    │
┌───────────────────────────────────▼────────────────────────────────────┐
│                       SHARED CORE INFRASTRUCTURE                       │
│  shared.io      shared.logging     shared.visualization  shared.report │
└───────────────────────────────────┬────────────────────────────────────┘
                                    │
┌───────────────────────────────────▼────────────────────────────────────┐
│                    SCIENTIFIC ENGINES & DATA STORES                    │
│   Biopython   •   NumPy / SciPy / Statsmodels   •   RCSB / NCBI / GEO  │
└────────────────────────────────────────────────────────────────────────┘
```

---

## 2. Directory Layout & Module Responsibilities

```
bioinformatics-portfolio/
├── shared/
│   ├── io/                 # FASTA, PDB, PDBQT, JSON, CSV parsers and writers
│   ├── logging/            # Structured logging with Rich and rotating handlers
│   ├── visualization/      # Publication styling (300 DPI, colorblind-safe palettes)
│   ├── reporting/          # Standalone Jinja2 HTML report generator with inline CSS/SVG
│   └── utils/              # API rate-limiters, hash verification, biological constants
│
└── projects/
    ├── 01-sequence-homology/
    │   ├── src/            # bioseq CLI, sequence statistics, BLAST engine, report builder
    │   ├── tests/          # Pytest suite with mock BLAST fixtures
    │   ├── examples/       # Verified FASTA files (Androgen Receptor, Insulin)
    │   ├── configs/        # Configurable BLAST and filtering thresholds
    │   └── results/        # Generated JSON, CSV, SVG, and HTML outputs
    │
    ├── 02-androgen-receptor-phylogenetics/
    │   ├── src/            # Alignment parser, entropy calculator, tree reconstructor
    │   ├── data/           # Curated FASTA dataset across 12 vertebrate orthologs
    │   ├── workflows/      # Snakemake DAG and parameter file
    │   ├── figures/        # Newick tree plots, conservation heatmaps
    │   └── results/        # Alignments (.aln, .fasta), Newick (.nwk), summary JSON
    │
    ├── 03-alphafold-structure-analysis/
    │   ├── src/            # PDB/mmCIF parser, pLDDT/PAE calculator, DSSP/SASA engine
    │   ├── data/           # Human AR LBD PDB (1E3G) and AlphaFold model (P10275)
    │   ├── scripts/        # Automated PyMOL (.pml) and ChimeraX (.cxc) sessions
    │   ├── figures/        # pLDDT residue profile, PAE heatmap, pocket cross-sections
    │   └── results/        # Structural assessment tables & summary report
    │
    ├── 04-protein-ligand-docking/
    │   ├── src/            # Receptor/ligand prep, Vina wrapper, PLIP contact analyzer
    │   ├── data/           # Cleaned PDB targets, 3D SDF ligand structures (DHT, Enz)
    │   ├── configs/        # Search grid coordinates, exhaustiveness, seeds
    │   ├── workflows/      # Dockerfile & Snakemake pipeline
    │   ├── figures/        # 2D/3D interaction diagrams, binding energy bar charts
    │   └── results/        # Posed PDBQTs, binding affinities (kcal/mol), contact CSV
    │
    └── 05-rna-seq/
        ├── src/            # Normalization, Negative Binomial DE, PCA, ORA pathway engine
        ├── data/           # Real GEO count matrix & sample metadata (GSE153664)
        ├── configs/        # Experimental contrasts, FDR thresholds, design formulas
        ├── workflows/      # Snakemake orchestration pipeline
        ├── figures/        # PCA plot, Volcano plot, MA plot, expression heatmap, GO bars
        └── results/        # DESeq2-equivalent DE tables, top gene tables, HTML report
```

---

## 3. Engineering Practices

- **Strict Type Hints**: Full static type coverage enforced with `mypy --disallow-untyped-defs`.
- **Fault-Tolerant Network Calls**: Remote APIs (NCBI, UniProt, RCSB PDB) implement retry logic with exponential jitter.
- **Self-Contained Standalone HTML Reports**: Reports include embedded SVG graphics, modern CSS layouts, and tabular data tables without external runtime CDN dependencies.
- **Fail-Safe Fallbacks**: Tree reconstruction and structure parsing include pure-Python native fallbacks if external command-line binaries (e.g., IQ-TREE or DSSP) are not present in the environment.
