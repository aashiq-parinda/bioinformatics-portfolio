# 🧬 Bioinformatics & Computational Biology Freelance Services & Rate Card

**Ashraf Khan** | *Senior Bioinformatics Software Engineer & Systems Architect*  
📍 Mumbai, India | 🌐 [Portfolio](https://github.com/aashiq-parinda/bioinformatics-portfolio) | ✉️ [ashrafk.salim1@gmail.com](mailto:ashrafk.salim1@gmail.com) | 📱 +91-8779559898

---

## 💡 Executive Value Proposition

Unlike academic bioinformaticians who deliver disconnected scripts or raw notebooks that break on other machines, I bring **5+ years of enterprise software engineering, HIPAA-compliant cloud architecture, and production-grade computational biology**.

Every deliverable is:
- ✅ **Fully Reproducible:** Containerized with Docker / Singularity, pinned dependencies, and deterministic random seeds.
- ✅ **Audit-Ready:** Accompanied by automated unit test suites (`pytest`), linting (`ruff`), and continuous integration.
- ✅ **Publication & Board-Ready:** High-resolution 300+ DPI vector graphics, interactive HTML portals, and executive summaries.
- ✅ **Confidential & Compliant:** Strict NDA execution, isolated compute environments, and HIPAA/PHI compliance expertise.

---

## 💰 Engagement & Pricing Models Summary

| Engagement Model | US & Global Rate (USD) | Indian Domestic Rate (INR) | Best Suited For |
| :--- | :--- | :--- | :--- |
| **Hourly Consulting / Code Review** | **$85 – $120 / hr** | **₹3,500 – ₹5,000 / hr** | Ad-hoc debugging, pipeline troubleshooting, methodology review |
| **Part-Time Retainer (15 hrs/week)** | **$4,500 – $6,000 / month** | **₹1,80,000 – ₹2,50,000 / month** | Ongoing lab support, biotech startup fractional bioinformatician |
| **Fixed-Scope Milestone Project** | **$1,500 – $8,500 / project** | **₹60,000 – ₹3,50,000 / project** | Defined end-to-end deliverables with fixed turnaround |

---

## 📦 Detailed Service Catalog & Fixed-Price Tiers

### 1. Transcriptomics & Single-Cell Omics

#### Service 1.1: Bulk RNA-Seq Differential Expression & Pathway Analysis
- **Scope:** Raw FASTQ QC or count matrix processing, DESeq2 median-of-ratios normalization, dispersion estimation, Wald test / GLM modeling, Benjamini-Hochberg FDR correction.
- **Deliverables:**
  - Normalized count tables and filtered DE gene lists ($|\log_2\text{FC}| \ge 1$, $\text{FDR} < 0.05$).
  - Publication-ready PCA scatter, Volcano plots, clustered expression heatmaps.
  - Over-Representation Analysis (ORA) & GSEA across GO terms and KEGG pathways.
  - Standalone interactive HTML report.
- **Turnaround:** 3–5 Business Days.
- **Pricing:** **$1,500 USD** | **₹60,000 INR** *(up to 12 samples)*

#### Service 1.2: Single-Cell RNA-Seq (scRNA-seq) End-to-End Analysis
- **Scope:** 10x Chromium count matrix ingestion, QC filtering (UMI, mito/ribo %, doublets via `Scrublet`), normalization, HVG selection, batch correction (`Harmony`/`scVI`), Leiden clustering, and marker gene annotation.
- **Deliverables:**
  - Annotated `AnnData` (`.h5ad`) object ready for exploration.
  - 2D UMAP visualizations, marker gene dotplots, cell-type composition bar charts.
  - Cluster differential expression tables.
- **Turnaround:** 5–7 Business Days.
- **Pricing:** **$2,800 USD** | **₹1,15,000 INR** *(up to 40,000 cells / 4 samples)*

#### Service 1.3: Spatial Transcriptomics & Tumor Microenvironment Deconvolution
- **Scope:** 10x Visium spatial slide alignment with matched scRNA-seq reference; spot-level quality control; Bayesian probabilistic cell-type deconvolution (`Cell2location`); spatial autocorrelation (Moran's $I$ SVGs); paracrine ligand-receptor cell-cell communication.
- **Deliverables:**
  - Deconvolved cell density overlays mapped onto high-resolution H&E images.
  - Spatial interaction network graphs and prioritized ligand-receptor pairs.
  - Complete reproducible Python/Squidpy workflow script.
- **Turnaround:** 7–10 Business Days.
- **Pricing:** **$4,200 USD** | **₹1,75,000 INR**

---

### 2. AI for Science & Protein Engineering

#### Service 2.1: AlphaFold Structure Assessment & Binding Site Profiling
- **Scope:** Structural analysis of experimental PDB or AlphaFold predicted models; per-residue pLDDT and Predicted Aligned Error (PAE) domain boundary extraction; pocket cleft detection and solvent-accessible surface area (SASA) calculation.
- **Deliverables:**
  - High-resolution confidence metrics plots and residue tolerance profiles.
  - Automated PyMOL / ChimeraX ray-tracing session scripts.
  - Structural summary report.
- **Turnaround:** 2–3 Business Days.
- **Pricing:** **$1,200 USD** | **₹45,000 INR** *(per target protein)*

#### Service 2.2: ESM-2 Deep Mutational Scanning & Zero-Shot Variant Effect Prediction
- **Scope:** Full in-silico mutational landscape scoring across all $19 \times L$ single amino-acid substitutions using Meta's `ESM-2` masked language model; log-likelihood ratio (LLR) scoring; validation against experimental DMS benchmarks (MaveDB); identification of stabilizing vs. deleterious mutations.
- **Deliverables:**
  - Complete $19 \times L$ mutational fitness CSV matrix.
  - Publication-quality interactive 2D mutational heatmaps.
  - PDB file with B-factors mapped to mutational tolerance for 3D visualization.
- **Turnaround:** 3–5 Business Days.
- **Pricing:** **$2,500 USD** | **₹1,00,000 INR**

#### Service 2.3: Interactive 3D Protein Analytics Web Application (Streamlit / Gradio)
- **Scope:** Building a bespoke, cloud-deployable web application for your team or clients to upload FASTA/PDB files, run in-silico predictions, explore interactive heatmaps, and rotate 3D Mol* structures in real time.
- **Deliverables:**
  - Standalone Dockerized Streamlit/Gradio web application with clean UI/UX.
  - Automated deployment to AWS ECS, GCP Cloud Run, or Hugging Face Spaces.
- **Turnaround:** 7–12 Business Days.
- **Pricing:** **$3,800 USD** | **₹1,50,000 INR**

---

### 3. Computational Drug Discovery & Biophysics

#### Service 3.1: Virtual Screening & Molecular Docking Campaign
- **Scope:** Target receptor active-site preparation; ligand 3D conformer generation, protonation, and energy minimization (`RDKit`/`OpenBabel`); parameterized `AutoDock Vina` screening; binding affinity scoring ($\Delta G$ kcal/mol); non-covalent protein-ligand interaction fingerprinting (`PLIP`).
- **Deliverables:**
  - Ranked compound hits spreadsheet with binding energies and efficiency indices.
  - 2D/3D interaction diagrams (H-bonds, salt bridges, $\pi$-stacking).
  - PDBQT/PDB docked complexes ready for visual inspection.
- **Turnaround:** 4–6 Business Days.
- **Pricing:** **$2,200 USD** | **₹90,000 INR** *(up to 1,000 compound conformations)*

#### Service 3.2: Explicit-Solvent Molecular Dynamics (OpenMM) & MM-PBSA Free Energy
- **Scope:** Converting docked hits into all-atom explicit-solvent systems; ligand parameterization (`OpenFF 2.1` / `GAFF2`); TIP3P water box with 0.15 M NaCl; 4-stage GPU simulation protocol in `OpenMM` (Minimization $\to$ NVT $\to$ NPT $\to$ 20–50 ns Production); trajectory analysis (`MDAnalysis`: RMSD, RMSF, H-bond persistence); binding free energy estimation ($\Delta G_{\text{bind}}$) via MM-PBSA.
- **Deliverables:**
  - Complete time-series RMSD/RMSF and H-bond kinetics graphs.
  - Decomposed thermodynamic free energy report ($\Delta E_{\text{MM}}$, $\Delta G_{\text{solv}}$, $\Delta G_{\text{bind}}$).
  - Aligned trajectory file (`.xtc` / `.dcd`) and PyMOL animation render script.
- **Turnaround:** 7–12 Business Days.
- **Pricing:** **$4,500 USD** | **₹1,85,000 INR** *(per complex, 50 ns simulation)*

---

### 4. Production Pipeline Engineering & DevOps

#### Service 4.1: Production Nextflow DSL2 Workflow Development (nf-core Compliant)
- **Scope:** Converting legacy shell scripts or monolithic Python scripts into modular Nextflow DSL2 pipelines. Includes process containers (Docker/Singularity), automatic memory retry handlers, configurable profiles (`-profile docker,awsbatch,slurm`), and automated `MultiQC` synthesis.
- **Deliverables:**
  - Complete Nextflow DSL2 repository with modular `modules/` and `subworkflows/`.
  - Continuous integration (GitHub Actions) testing on synthetic minimal data.
  - Comprehensive documentation and run-command guide.
- **Turnaround:** 7–14 Business Days.
- **Pricing:** **$4,000 – $6,500 USD** | **₹1,60,000 – ₹2,60,000 INR**

#### Service 4.2: Custom Bioinformatics CLI / Python Package Development
- **Scope:** Engineering an open-source or proprietary command-line tool (CLI) using `Typer` and `Rich`; automated argument parsing, robust error handling, fast I/O parsers (FASTA, VCF, PDB), logging, and complete `pytest` coverage.
- **Deliverables:**
  - Pip/Conda installable Python package (`pyproject.toml`).
  - Automated test suite with 85%+ test coverage.
  - Markdown / Sphinx documentation.
- **Turnaround:** 5–8 Business Days.
- **Pricing:** **$2,500 USD** | **₹1,00,000 INR**

---

### 5. Clinical Genomics & Precision Oncology

#### Service 5.1: Somatic Variant Curation & AMP/ASCO/CAP 4-Tier Engine
- **Scope:** Ingestion and normalization of tumor-normal somatic VCFs; functional consequence annotation (`Ensembl VEP`); population frequency filtering (`gnomAD`); automated querying of clinical oncology knowledgebases (`CIViC API`, `ClinVar`); deterministic classification into **Tier I** (FDA-approved/standard of care) through **Tier IV** (benign).
- **Deliverables:**
  - Structured JSON and CSV tables of tiered somatic variants.
  - Actionable drug and clinical trial matching matrix.
  - Automated physician-ready molecular pathology PDF report.
- **Turnaround:** 5–8 Business Days.
- **Pricing:** **$3,500 USD** | **₹1,40,000 INR**

---

## 🛡️ Terms of Service, Security & Governance

1. **Non-Disclosure Agreement (NDA):** Signed prior to receiving any raw sequencing, chemical structures, or patient-derived data.
2. **Data Custody & Security:** All computation is performed on encrypted, ephemeral instances (AWS / GCP) adhering to HIPAA guidelines. No client data is ever uploaded to unauthorized third-party public AI models.
3. **Intellectual Property (IP):** 100% of code, models, analyses, and deliverables belong to the client upon milestone payment completion.
4. **Payment Schedule:**
   - Fixed-price projects: 50% deposit upon contract initiation, 50% upon delivery and client validation.
   - Retainers: Invoiced at the beginning of each 30-day billing cycle.
   - Payment rails supported: Wise (USD/EUR/GBP), Stripe, Direct Wire Transfer, and Razorpay/UPI/NEFT (INR).

---

## 📬 How to Book a Project

1. **Send an inquiry:** Email `ashrafk.salim1@gmail.com` with your project scope, dataset size, and timeline requirements.
2. **Introductory Scoping Call (30 mins):** We review biological hypotheses, technical feasibility, compute environment, and agree on deliverables.
3. **Formal Proposal & Contract:** Detailed Statement of Work (SOW) with fixed milestones, deliverables, and dates provided within 24 hours.
