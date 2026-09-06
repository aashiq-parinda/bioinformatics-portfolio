# Scientific Tooling & Algorithmic Inventory

This inventory documents the scientific software, mathematical formulations, and algorithmic foundations implemented across the repository.

---

## 1. Sequence & Homology Tooling (`projects/01-sequence-homology`)

- **Biopython (`Bio.SeqUtils.ProtParam`)**:
  - Isoelectric Point Calculation: Bisection algorithm over pKa values (Bjellqvist / EMBOSS scales).
  - Molecular Weight: Monoisotopic / average isotopic mass summation minus water condensation.
  - GRAVY (Grand Average of Hydropathy): Kyte-Doolittle hydropathy index mean $\frac{1}{N}\sum_{i=1}^N H_i$.
  - Instability Index: Guruprasad et al. (1990) dipeptide instability weight matrix.
- **BLAST+ & NCBI QBLAST API**:
  - Execution of local `blastp` binary and remote `NCBIWWW.qblast` with automated exponential backoff and XML/JSON hit extraction.
  - Scoring matrices: BLOSUM62, BLOSUM45, PAM250, PAM70 with affine gap penalties.

---

## 2. Phylogenetic Inference & Multiple Alignment (`projects/02-androgen-receptor-phylogenetics`)

- **Alignment Algorithms**:
  - **MAFFT** (L-INS-i / FFT-NS-2 algorithms) & **MUSCLE** (iterative profile-profile alignment).
  - Shannon Entropy Conservation Scoring:
    $$H_i = -\sum_{a \in \mathcal{A}} P(a) \log_2 P(a)$$
- **Phylogenetic Reconstruction**:
  - **IQ-TREE / FastTree / Phylo (Neighbor-Joining & UPGMA)** with JTT/WAG substitution models and $\Gamma$-distributed rate heterogeneity across sites.
  - Felsenstein Non-Parametric Bootstrap (1,000 replicates) for branch support estimation.

---

## 3. Structural Biology & AlphaFold Evaluation (`projects/03-alphafold-structure-analysis`)

- **AlphaFold v2 Metrics**:
  - **pLDDT (Predicted Local Distance Difference Test)**: Per-residue confidence scores ($0-100$) categorized into Very High ($>90$), Confident ($70-90$), Low ($50-70$), and Very Low ($<50$, indicative of intrinsic disorder).
  - **PAE (Predicted Aligned Error)**: Inter-domain rigid-body distance error matrix ($\text{\AA}$).
- **Structure Analytics**:
  - Secondary Structure Assignment: DSSP hydrogen bonding pattern estimation ($\alpha$-helices, $\beta$-sheets, coils).
  - SASA (Solvent Accessible Surface Area): Shrake-Rupley / Lee-Richards rolling sphere algorithm ($r_{\text{probe}} = 1.4\text{ \AA}$).
  - PyMOL & ChimeraX automation scripts for ray-traced rendering of active sites and ligand binding pockets.

---

## 4. Molecular Docking & Interaction Profiling (`projects/04-protein-ligand-docking`)

- **AutoDock Vina Engine**:
  - Iterated Local Search with BFGS quasi-Newton local optimization.
  - Empirical free energy scoring function incorporating steric, hydrophobic, hydrogen bonding, and torsional penalty terms:
    $$\Delta G_{\text{bind}} = \sum \text{vdw} + \text{hbond} + \text{elec} + \text{desolv} + \text{tors}$$
- **RDKit Cheminformatics**:
  - Ligand 3D conformer generation (ETKDGv3 method) and MMFF94 force-field energy minimization.
  - Gasteiger-Marsili partial charge assignment.
- **PLIP Interaction Profiling**:
  - Identification of hydrogen bonds (geometric angle $> 100^\circ$, distance $< 3.5\text{ \AA}$), salt bridges, $\pi$-stacking, and hydrophobic contact shells.

---

## 5. Transcriptomics & Differential Expression (`projects/05-rna-seq`)

- **Statistical Modeling**:
  - Negative Binomial generalized linear model:
    $$K_{ij} \sim \text{NB}(\mu_{ij}, \alpha_i), \quad \mu_{ij} = s_j q_{ij}$$
  - Size factor estimation via Median-of-Ratios normalization.
  - Dispersion estimation via empirical Bayes shrinkage toward a trended curve.
  - Wald hypothesis testing with Benjamini-Hochberg False Discovery Rate (FDR / adjusted $p$-value) correction.
- **Functional Enrichment (ORA)**:
  - Hypergeometric distribution test (Fisher's exact test) against Gene Ontology (GO Biological Process, Molecular Function, Cellular Component) and KEGG pathways.
