# Professional Bioinformatics Freelancing & Client Services Matrix

This guide outlines professional bioinformatics services mapped directly to the competencies demonstrated across this portfolio. Each service tier defines client deliverables, typical scopes of work, and honest boundaries distinguishing beginner-ready engagements from advanced senior-only contracts.

---

## 1. Direct Service Mapping by Portfolio Flagship

```
┌────────────────────────────────────────────────────────────────────────┐
│               COMPUTATIONAL BIOLOGY CLIENT SERVICES                    │
├───────────────────┬────────────────────────────────────────────────────┤
│ 01. Sequence CLI  │ Automated Bio-Pipelines • Homology • Batch BLAST   │
│ 02. Phylogenetics │ Ortholog Screening • Multiple Alignment • Tree Evo │
│ 03. AlphaFold 3D  │ Structure Assessment • pLDDT/PAE • Binding Pockets │
│ 04. Docking Vina  │ Virtual Screening • Ligand Posing • Contact Profil │
│ 05. RNA-Seq DE    │ End-to-End Counts to Volcano • Pathway Enrichment  │
└───────────────────┴────────────────────────────────────────────────────┘
```

---

## 2. Core Service Offerings (Immediate Client Readiness)

### Service 1: Custom Bioinformatics CLI & Automation Pipelines
*Demonstrated by Project 01 (`bioseq`)*
- **Client Need**: Laboratories and biotech startups burdened with repetitive manual NCBI searches, sequence formatting, and ad-hoc Excel manipulations.
- **Deliverables**:
  - Python/Typer CLI installable via pip/conda.
  - Robust batch FASTA validation, translation, and physicochemical profiling.
  - Automated local BLAST+ orchestration and filtering into JSON/CSV and executive HTML summaries.
- **Scope**: Sequence automation, data wrangling, API integration, and parsing utilities.

### Service 2: Multiple Sequence Alignment & Evolutionary Phylogenetics
*Demonstrated by Project 02 (`androgen-receptor-phylogenetics`)*
- **Client Need**: Evolutionary biology researchers and antibody engineering teams requiring ortholog profiling, active-site conservation analysis, and publication-ready phylogenetic trees.
- **Deliverables**:
  - Curated ortholog sequence retrieval from UniProt/NCBI.
  - MAFFT/MUSCLE alignment with position-specific Shannon entropy calculations.
  - Maximum-Likelihood / Neighbor-Joining phylogenetic tree reconstruction with bootstrap validation.
  - Vector graphics (SVG/PDF) for publication submission.

### Service 3: Protein Structural Characterization & AlphaFold Analysis
*Demonstrated by Project 03 (`alphafold-structure-analysis`)*
- **Client Need**: Drug discovery teams needing independent validation of predicted structures, domain boundary determination, and secondary structure quantification.
- **Deliverables**:
  - AlphaFold v2 / v3 confidence profiling (pLDDT per residue and PAE domain matrices).
  - Binding cleft and active-site residue mapping.
  - Automated PyMOL and ChimeraX visualization scripts for presentation to wet-lab stakeholders.

### Service 4: Reproducible Molecular Docking & Interaction Profiling
*Demonstrated by Project 04 (`protein-ligand-docking`)*
- **Client Need**: Small-molecule discovery teams requiring initial computational screening of compound libraries against verified target structures.
- **Deliverables**:
  - High-throughput receptor and ligand preparation (RDKit / PDBQT generation).
  - Parameterized AutoDock Vina execution with search space optimization.
  - Pose clustering, binding affinity ranking ($\Delta G$ kcal/mol), and PLIP interaction fingerprints.
  - Scientific report with explicit computational disclaimers.

### Service 5: End-to-End RNA-Seq Differential Expression & Pathway Analysis
*Demonstrated by Project 05 (`rna-seq`)*
- **Client Need**: Research labs with raw or count-level RNA-seq data needing statistically sound differential expression and biological interpretation.
- **Deliverables**:
  - Quality assessment (FastQC/MultiQC metric synthesis).
  - Count normalization (DESeq2 median-of-ratios methodology).
  - Statistical testing (Wald test, Benjamini-Hochberg FDR).
  - Interactive/Publication figures: Sample PCA, Volcano plots, MA plots, Hierarchical heatmaps.
  - Functional enrichment: GO terms (BP, MF, CC) and KEGG pathway over-representation analysis.

---

## 3. Advanced Services Requiring Extended Domain Expertise

To maintain rigorous scientific ethics and client trust, the following services are explicitly categorized as **Advanced** and should only be undertaken in collaboration with specialized domain experts or after extensive wet-lab validation experience:

| Advanced Service | Technical Requirements Beyond Standard Pipelines |
| :--- | :--- |
| **All-Atom Molecular Dynamics (MD)** | Microsecond Gromacs/AMBER simulation, force-field parameterization for novel ligands, free-energy perturbation (FEP). |
| **Clinical Diagnostic Variant Curation** | ACMG/AMP clinical variant guidelines, clinical trial compliance, FDA/EMA regulatory reporting. |
| **De Novo Protein Design (RFdiffusion / ProteinMPNN)** | Deep generative modeling, experimental yeast/phage display validation protocols. |
| **Single-Cell Trajectory & Spatial Multi-Omics** | Velocity modeling (scVelo), spatial transcriptomics deconvolution, batch integration across heterogeneous sequencing platforms. |

---

## 4. Ethical & Client Communication Standards

1. **Clear Distinction of Hypotheses**: Computational docking and virtual screening results are delivered as prioritized hypotheses for wet-lab assay validation, never as clinical efficacy claims.
2. **Deterministic & Audit-Ready Deliverables**: All code handed over to clients includes unit tests, environment lockfiles, reproducible seeds, and clear README execution steps.
3. **Transparent Data Governance**: Client proprietary sequence data is handled in isolated secure environments and never transmitted over unauthorized public APIs.
