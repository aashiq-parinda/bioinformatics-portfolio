# Project 08: Production Pathogen Genomic Surveillance Pipeline (Nextflow DSL2)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-brightgreen.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/container-Docker%20%7C%20Singularity-blue.svg)](https://www.docker.com/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![MultiQC](https://img.shields.io/badge/reporting-MultiQC-orange.svg)](https://multiqc.info/)

## Overview

Modern infectious disease surveillance, outbreak tracking, and antimicrobial resistance (AMR) monitoring require robust, automated, and horizontally scalable genomics pipelines. Ad-hoc shell scripts fail across heterogeneous compute infrastructure and cannot provide the auditability required in clinical and public-health laboratories.

This project delivers an enterprise-grade **Pathogen Genomic Surveillance Pipeline** architected strictly using **Nextflow DSL2** adhering to **nf-core** development standards. The pipeline ingests Illumina paired-end short reads or Oxford Nanopore (ONT) long reads, executes automated quality control, reference alignment, high-accuracy variant calling, consensus genome generation, clade/lineage typing, and AMR gene profiling, culminating in an interactive **MultiQC** synthesis report.

---

## Architecture & Workflow Graph

```
                   INPUT READS (Illumina / Oxford Nanopore)
                                     │
                     ┌───────────────┴───────────────┐
                     ▼                               ▼
             QC & TRIMMING                     QC & FILTERING
        (FastQC + fastp [Illumina])         (NanoPlot + NanoFilt [ONT])
                     │                               │
                     └───────────────┬───────────────┘
                                     ▼
                            REFERENCE ALIGNMENT
                     (BWA-MEM2 / Minimap2 + Samtools)
                                     │
                     ┌───────────────┴───────────────┐
                     ▼                               ▼
              VARIANT CALLING                 DE NOVO / AMR
         (FreeBayes / Clair3 + bcftools)     (AMRFinderPlus / CARD RGI)
                     │                               │
                     ▼                               ▼
            CONSENSUS ASSEMBLY                AMR PROFILE &
        (Masked low-coverage Ns)             RESISTANCE CALLS
                     │                               │
                     ▼                               │
             LINEAGE TYPING                          │
         (Nextclade / Pangolin / MLST)               │
                     │                               │
                     └───────────────┬───────────────┘
                                     ▼
                          MULTIQC SYNTHESIS REPORT
                     (Unified interactive HTML dashboard)
```

---

## Key Features

1. **Modular DSL2 Architecture (`modules/` & `subworkflows/`):**
   - Encapsulated, reusable Nextflow DSL2 process modules with version-tagged container directives.
   - Decoupled subworkflows for read preprocessing, alignment, variant calling, and pathogen profiling.

2. **Dual Sequencing Technology Support:**
   - **Illumina:** Ultra-fast quality trimming via `fastp`, indexing and alignment via `bwa-mem2`, and diploid/haploid variant calling via `freebayes`/`bcftools`.
   - **Oxford Nanopore (ONT):** Length/quality filtering via `NanoFilt`, long-read alignment via `minimap2`, and deep-learning variant calling via `clair3`.

3. **Antimicrobial Resistance & Clade Profiling (`modules/amr.nf`):**
   - Automated screening for acquired resistance genes, point mutations, and virulence determinants via NCBI `AMRFinderPlus` and the Comprehensive Antibiotic Resistance Database (`CARD RGI`).
   - Lineage assignment for viral surveillance (e.g. Nextclade / Pangolin) and bacterial Multi-Locus Sequence Typing (`MLST`).

4. **Production Portability & Cloud Execution:**
   - Multi-environment profiles in `nextflow.config`: `-profile docker`, `-profile singularity`, `-profile conda`, `-profile slurm`, and `-profile awsbatch`.
   - Automated memory/CPU retry handlers on out-of-memory errors (`errorStrategy = 'retry'`).

5. **Automated Quality Synthesis:**
   - Aggregates trimming metrics, alignment statistics, variant counts, and AMR tables into a single publication-ready `MultiQC` HTML report.

---

## Directory Structure

```
08-nextflow-pathogen-surveillance/
├── assets/
│   ├── multiqc_config.yaml       # MultiQC branding & custom metrics tables
│   └── test_samplesheet.csv      # Test sample manifest (sample, fastq_1, fastq_2, platform)
├── configs/
│   ├── base.config               # Compute resource allocations & retry rules
│   ├── docker.config             # Pinned container images per tool
│   └── test.config               # Minimal configuration for quick smoke testing
├── data/
│   ├── references/               # Target reference genomes (e.g., SARS-CoV-2, M. tuberculosis)
│   └── raw_fastq/                # Input FASTQ reads
├── docs/
│   ├── parameters.md             # CLI parameter specification
│   └── output_schema.md          # Description of all emitted files
├── modules/
│   ├── alignment.nf              # BWA-MEM2 / Minimap2 alignment modules
│   ├── amr_profiling.nf          # AMRFinderPlus / CARD modules
│   ├── qc_trimming.nf            # FastQC, NanoPlot, and fastp modules
│   ├── report.nf                 # MultiQC module
│   └── variant_calling.nf        # FreeBayes, bcftools, consensus modules
├── subworkflows/
│   ├── align_and_call.nf         # Combined alignment and variant detection
│   └── profile_pathogen.nf       # Clade typing and resistance profiling
├── nextflow.config               # Master Nextflow configuration
├── main.nf                       # Pipeline entrypoint
└── workflows/
    └── pathogen_surveillance.nf  # Primary workflow logic
```

---

## Quick Start

```bash
# 1. Run pipeline on test dataset using Docker
nextflow run projects/08-nextflow-pathogen-surveillance/main.nf \
    -profile docker,test \
    --input projects/08-nextflow-pathogen-surveillance/assets/test_samplesheet.csv \
    --outdir projects/08-nextflow-pathogen-surveillance/results/

# 2. View generated MultiQC Report
open projects/08-nextflow-pathogen-surveillance/results/multiqc/multiqc_report.html
```
