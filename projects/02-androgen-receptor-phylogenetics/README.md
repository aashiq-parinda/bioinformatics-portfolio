# Project 02: Androgen-Receptor Comparative Phylogenetic Analysis

[![Workflow: Snakemake](https://img.shields.io/badge/Workflow-Snakemake-blue.svg)](workflows/Snakefile)
[![Tests: Pytest](https://img.shields.io/badge/tests-passing-brightgreen.svg)](tests/)

A reproducible comparative evolutionary genomics pipeline analyzing androgen receptor (*AR* / *NR3C4*) sequence conservation and divergence across 10 representative vertebrate species.

---

## 🧬 Biological & Evolutionary Context

Nuclear steroid receptors evolved through ancestral whole-genome duplications. The androgen receptor coordinates transcriptional programs via:
1. **DNA-Binding Domain (DBD)**: Two zinc-finger motifs recognizing androgen response elements (AREs: `5'-AGAACA-3'`).
2. **Hinge Region**: Flexible linker containing nuclear localization signals.
3. **Ligand-Binding Domain (LBD)**: A 12-$\alpha$-helix sandwich forming the hydrophobic steroid hormone pocket.

This project reconstructs vertebrate receptor divergence and quantifies domain-specific evolutionary constraints using **Shannon Entropy**:
$$H_i = -\sum_{a \in \mathcal{A}} P(a) \log_2 P(a)$$

> **Important Domain Safety Notice:** This is an evolutionary and computational biology analysis. It contains NO steroid cycle recommendations, anabolic dosing instructions, or clinical treatment advice.

---

## 🔬 Taxonomic Dataset

The analysis evaluates verified UniProt/NCBI records spanning major vertebrate classes:
- **Mammalia**: Human (`P10275`), Chimpanzee (`Q95180`), Rhesus macaque (`P51124`), Dog (`Q95197`), Bovine (`Q9TT90`), Mouse (`P19091`), Rat (`P15207`).
- **Aves**: Chicken (`P79780`).
- **Amphibia**: African clawed frog (`Q91843`).
- **Actinopterygii**: Zebrafish (`Q90WF7`, basal outgroup).

---

## 🛠️ Execution Instructions

### Run via Python Pipeline
```bash
python3 projects/02-androgen-receptor-phylogenetics/src/run_pipeline.py \
  --config projects/02-androgen-receptor-phylogenetics/configs/pipeline_config.json
```

### Run via Snakemake Workflow
```bash
snakemake --cores 2 -s projects/02-androgen-receptor-phylogenetics/workflows/Snakefile
```

### Run Automated Tests
```bash
pytest projects/02-androgen-receptor-phylogenetics/tests/ -v
```

---

## 📊 Deliverables & Generated Artifacts

- `results/ar_vertebrates_aligned.fasta`: Clean multiple sequence alignment.
- `results/per_residue_conservation.csv`: Column-by-column Shannon entropy and consensus frequencies.
- `results/ar_phylogenetic_tree.nwk`: Valid Newick phylogenetic tree format.
- `results/phylogenetic_analysis_report.html`: Standalone interactive HTML report with species cards and stats.
- `figures/phylogenetic_tree.png`: Publication-grade phylogram with taxon labels.
- `figures/conservation_profile.png`: Per-residue entropy plot with functional domain overlays.
- `figures/domain_conservation.png`: Comparison of mean evolutionary entropy across domains.
- `notebooks/02_phylogenetics_walkthrough.ipynb`: Interactive Jupyter walkthrough notebook.
