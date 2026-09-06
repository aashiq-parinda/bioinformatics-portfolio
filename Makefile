.PHONY: help install install-dev test lint format check clean run-all

PYTHON ?= python3
PIP ?= $(PYTHON) -m pip
PYTEST ?= $(PYTHON) -m pytest
RUFF ?= $(PYTHON) -m ruff

help:
	@echo "Available commands:"
	@echo "  install        Install runtime dependencies"
	@echo "  install-dev    Install runtime and developer dependencies"
	@echo "  test           Run all pytest test suites"
	@echo "  lint           Run Ruff linter and static checks"
	@echo "  format         Run Ruff formatter"
	@echo "  check          Run formatting check, linting, and tests"
	@echo "  clean          Remove build and test cache artifacts"
	@echo "  run-p01        Run Project 01 (Sequence & Homology Analysis)"
	@echo "  run-p02        Run Project 02 (Androgen Receptor Phylogenetics)"
	@echo "  run-p03        Run Project 03 (AlphaFold Structure Analysis)"
	@echo "  run-p04        Run Project 04 (Protein-Ligand Docking Pipeline)"
	@echo "  run-p05        Run Project 05 (End-to-End RNA-Seq Pipeline)"
	@echo "  run-all        Run all 5 flagship project pipelines sequentially"

install:
	$(PIP) install -e .

install-dev:
	$(PIP) install -e ".[dev,workflow,docking]"

test:
	$(PYTEST) -v

lint:
	$(RUFF) check .

format:
	$(RUFF) format .

check: lint test

clean:
	rm -rf build/ dist/ *.egg-info .pytest_cache .coverage htmlcov/ .ruff_cache
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

run-p01:
	$(PYTHON) projects/01-sequence-homology/src/cli.py analyze projects/01-sequence-homology/examples/human_androgen_receptor.fasta --output projects/01-sequence-homology/results/

run-p02:
	$(PYTHON) projects/02-androgen-receptor-phylogenetics/src/run_pipeline.py --config projects/02-androgen-receptor-phylogenetics/configs/pipeline_config.json

run-p03:
	$(PYTHON) projects/03-alphafold-structure-analysis/src/run_analysis.py --config projects/03-alphafold-structure-analysis/configs/analysis_config.json

run-p04:
	$(PYTHON) projects/04-protein-ligand-docking/src/run_docking_pipeline.py --config projects/04-protein-ligand-docking/configs/docking_config.json

run-p05:
	$(PYTHON) projects/05-rna-seq/src/run_rnaseq_pipeline.py --config projects/05-rna-seq/configs/rnaseq_config.json

run-all: run-p01 run-p02 run-p03 run-p04 run-p05
