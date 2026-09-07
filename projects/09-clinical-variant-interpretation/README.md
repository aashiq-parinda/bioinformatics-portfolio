# Project 09: Clinical Somatic Variant Interpretation & Precision Oncology Engine

[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/)
[![Guidelines: AMP/ASCO/CAP](https://img.shields.io/badge/guidelines-AMP%2FASCO%2FCAP%20Tiers-red.svg)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5707196/)
[![API: CIViC](https://img.shields.io/badge/database-CIViC%20%7C%20ClinVar-blue.svg)](https://civicdb.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

In diagnostic precision oncology, calling somatic mutations from tumor-normal whole-exome or targeted next-generation sequencing (NGS) panels produces raw Variant Call Format (VCF) files with hundreds of candidate alterations. The principal clinical bottleneck is translating these raw genomic variants into actionable clinical intelligence: identifying oncogenic drivers, classifying therapeutic actionability according to formal medical guidelines, matching patients to approved drugs or open clinical trials, and generating an auditable, physician-ready molecular pathology report.

This project delivers an automated **Clinical Variant Interpretation & Drug Matching Engine** adhering to the consensus **AMP/ASCO/CAP (Association for Molecular Pathology / American Society of Clinical Oncology / College of American Pathologists)** four-tier somatic variant classification system.

---

## AMP/ASCO/CAP 4-Tier Somatic Classification System

```
┌───────────────────────────────────────────────────────────────────────────┐
│ TIER I: Variants of Strong Clinical Significance (Level A & B Evidence)   │
│ - FDA-approved therapies for specific tumor type (e.g. BRAF V600E in Mel) │
│ - Well-powered professional guidelines (NCCN, ASCO)                       │
├───────────────────────────────────────────────────────────────────────────┤
│ TIER II: Variants of Potential Clinical Significance (Level C & D)        │
│ - FDA-approved therapies in other tumor types (off-label)                 │
│ - Biomarkers for active Phase I/II/III clinical trial eligibility         │
├───────────────────────────────────────────────────────────────────────────┤
│ TIER III: Variants of Uncertain Significance (VUS)                        │
│ - Somatic alterations without convincing diagnostic/prognostic evidence   │
├───────────────────────────────────────────────────────────────────────────┤
│ TIER IV: Benign or Likely Benign Variants                                 │
│ - High population allele frequency in gnomAD (>1%) without cancer link    │
└───────────────────────────────────────────────────────────────────────────┘
```

---

## Key Features

1. **High-Performance VCF Normalization & Parsing (`src/vcf_parser.py`):**
   - Streamed parsing of multi-sample tumor-normal VCFs using `cyvcf2` and `pysam`.
   - Normalization: multi-allelic splitting, left-alignment of indels, and variant effect annotation filtering.
   - Quality filters: read depth ($DP \ge 50$), variant allele fraction ($VAF \ge 0.05$), strand bias, and population filtering via `gnomAD` ($AF < 0.01$ to eliminate germline polymorphisms).

2. **Functional & Pathogenicity Annotation (`src/annotator.py`):**
   - Variant effect integration using Ensembl VEP (transcript impact, amino acid alteration, consequence type).
   - In-silico damage prediction scores (CADD Phred, REVEL, SIFT, PolyPhen-2).
   - Clinical registry queries across NCBI `ClinVar` and `COSMIC` somatic recurrence tallies.

3. **Knowledgebase Integration & Evidence Grading (`src/civic_matcher.py`):**
   - Automated REST API integration with **CIViC** (Clinical Interpretation of Variants in Cancer).
   - Querying molecular evidence items (EID), clinical significance (Sensitivity/Response, Resistance, Adverse Response, Diagnostic, Prognostic), and disease ontology mapping.
   - Automated tier calculation mapping variants into **Tier I**, **Tier II**, **Tier III (VUS)**, or **Tier IV (Benign)** based on evidence strength and tumor histology.

4. **Targeted Drug & Trial Matching Engine (`src/drug_matcher.py`):**
   - Direct matching to FDA-approved precision therapeutics (e.g., Osimertinib for *EGFR* T790M, Sotorasib/Adagrasib for *KRAS* G12C, Trastuzumab for *ERBB2* amplification, Olaparib for *BRCA1/2* deleterious mutations).
   - Resistance mutation flagging (e.g., *EGFR* C797S resistance to 3rd-generation TKIs, *AR* L702H resistance to enzalutamide).

5. **Automated Clinical Molecular Pathology Report (`src/report_generator.py`):**
   - Generates an executive, audit-ready clinical PDF report and interactive HTML portal.
   - Includes patient metadata, tumor type, sequencing metrics, tiered mutation summary table, actionable drug indications, and referenced PubMed literature citations.

---

## Directory Structure

```
09-clinical-variant-interpretation/
├── configs/
│   └── clinical_config.yaml      # Quality thresholds, VAF cutoffs & API credentials
├── data/
│   ├── reference/                # Bed intervals, gene panels, gnomAD frequency files
│   └── vcfs/                     # Example somatic cancer VCFs (NSCLC, Melanoma, Colorectal)
├── notebooks/
│   └── 09_clinical_interpretation_walkthrough.ipynb # Interactive curation tutorial
├── reports/
│   └── sample_clinical_report.pdf # Generated sample physician PDF report
├── results/
│   ├── tiered_variants.json      # Structured JSON representation of tiered variants
│   ├── clinical_actionability.csv # Filtered actionable drug table
│   └── clinical_report.html      # Interactive clinical report
├── src/
│   ├── __init__.py
│   ├── vcf_parser.py             # cyvcf2 normalization, VAF calculation & QC
│   ├── annotator.py              # Functional impact & population allele frequency
│   ├── civic_matcher.py          # CIViC API evidence integration & tier assignment
│   ├── drug_matcher.py           # Precision oncology drug & trial matching
│   └── report_generator.py       # Automated PDF/HTML molecular tumor board report
└── tests/
    ├── test_vcf_parser.py        # Unit tests for VCF filtering & VAF thresholds
    ├── test_tier_assignment.py   # Validation of AMP/ASCO/CAP tier logic rules
    └── test_civic_api.py         # Mock and live tests for CIViC API endpoints
```

---

## Quick Start

```bash
# 1. Run somatic interpretation on patient tumor VCF
python3 projects/09-clinical-variant-interpretation/src/run_clinical_engine.py \
    --vcf projects/09-clinical-variant-interpretation/data/vcfs/patient_nsclc.vcf.gz \
    --tumor-type "Non-Small Cell Lung Carcinoma" \
    --output projects/09-clinical-variant-interpretation/results/

# 2. Run unit tests
pytest projects/09-clinical-variant-interpretation/tests/ -v

# 3. View the generated clinical report
open projects/09-clinical-variant-interpretation/results/clinical_report.html
```
