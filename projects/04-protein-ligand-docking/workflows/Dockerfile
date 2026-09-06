# Multi-stage reproducible Docker container for Molecular Docking
FROM continuumio/miniconda3:latest AS builder

WORKDIR /opt/docking

# Configure conda channels
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda config --set channel_priority strict

# Install system and scientific dependencies
RUN conda install -y \
    python=3.11 \
    autodock-vina \
    openbabel \
    rdkit \
    snakemake \
    biopython \
    pandas \
    scipy \
    matplotlib \
    seaborn \
    && conda clean -afy

# Copy repository
WORKDIR /app
COPY . /app

# Install package in editable mode
RUN pip install --no-deps -e .

ENV PYTHONUNBUFFERED=1
ENV RANDOM_SEED=42

# Default entrypoint runs Project 04 docking pipeline
CMD ["python3", "projects/04-protein-ligand-docking/src/run_docking_pipeline.py", "--config", "projects/04-protein-ligand-docking/configs/docking_config.json"]
