# Project 01: Automated Protein Sequence & Homology Analysis

[![Tool: bioseq](https://img.shields.io/badge/CLI-bioseq-blue.svg)](src/cli.py)
[![Test: Pytest](https://img.shields.io/badge/tests-passing-brightgreen.svg)](tests/)

A scientific CLI application and Python library for comprehensive physicochemical sequence profiling, automated local/remote BLAST homology querying, intelligent hit filtering, and automated HTML/JSON report generation.

---

## 🧬 Biological Context & Scope

Understanding primary protein sequences is the foundation of molecular biology. This pipeline calculates essential biophysical properties:
- **Isoelectric Point (pI)**: Theoretical pH at which the net charge equals zero (via Bjellqvist pKa scales).
- **Grand Average of Hydropathy (GRAVY)**: Kyte-Doolittle hydropathicity index indicating global membrane vs. soluble tendency.
- **Instability Index**: Guruprasad dipeptide statistical weights predicting in vivo protein stability.
- **Homology Search**: Automated BLAST+ / NCBI QBLAST querying with coverage/identity filtering to identify orthologs and paralogs.

---

## 🛠️ CLI Usage & Examples

### 1. Physicochemical Analysis
```bash
# Analyze Human Androgen Receptor (UniProt P10275)
bioseq analyze projects/01-sequence-homology/examples/human_androgen_receptor.fasta --output projects/01-sequence-homology/results/
```

### 2. Homology Search (NCBI QBLAST / Local BLAST+)
```bash
# Query NCBI SwissProt database with stringent E-value and identity thresholds
bioseq blast projects/01-sequence-homology/examples/human_insulin.fasta \
    --database swissprot \
    --evalue 1e-5 \
    --min-identity 30 \
    --min-coverage 50 \
    --output projects/01-sequence-homology/results/
```

---

## 📊 Pipeline Artifacts & Outputs

Each execution produces structured data and publication-quality figures:
- `<ID>_statistics.json`: Machine-readable metadata and numerical statistics.
- `<ID>_blast_hits.csv`: Tabular BLAST hits with query coverage, identity %, and bit scores.
- `<ID>_analysis_report.html`: Standalone responsive HTML summary report.
- `figures/amino_acid_composition.png`: High-DPI bar chart of residue percentages.
- `figures/hydropathy_profile.png`: Sliding-window Kyte-Doolittle hydropathy plot.
- `figures/blast_hits_distribution.png`: Scatter/bubble plot of identity vs. E-value.

---

## 🔬 Testing & Verification

```bash
# Run Project 01 pytest test suite
pytest projects/01-sequence-homology/tests/ -v
```
