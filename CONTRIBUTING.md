# Contributing Guidelines

Thank you for your interest in contributing to the **Bioinformatics Engineering Portfolio**.

This repository adheres to strict software engineering standards and biological data governance rules.

## Code Standards

1. **Python Version**: All code targets Python 3.11+.
2. **Type Hints**: All functions and public class methods must include complete PEP 484 type annotations.
3. **Linting & Formatting**: Enforced via [Ruff](https://astral.sh/ruff).
   ```bash
   make lint
   make format
   ```
4. **Testing**: All analytical logic, parsers, and statistical routines must be accompanied by `pytest` unit/integration tests with deterministic fixtures.
   ```bash
   make test
   ```
5. **Deterministic Execution**: Always fix random seeds (`RANDOM_SEED=42`) for any stochastic simulation or heuristic algorithm.

## Scientific Integrity & Safety Policy

- **No Medical/Cycle Advice**: Under no circumstances should scripts or documentation provide drug protocols, anabolic cycling advice, or clinical diagnoses.
- **Data Provenance**: Always cite public accession numbers (UniProt, PDB, NCBI GenBank, GEO/SRA) when introducing new biological examples.
- **Computational Disclaimers**: Computational docking results, homology models, and predictive algorithms must clearly declare their structural and statistical limitations.

## Pull Request Workflow

1. Fork the repository and create a descriptive feature branch (`feat/new-parser` or `fix/blast-timeout`).
2. Implement your changes with corresponding test cases.
3. Verify that `make check` passes cleanly.
4. Submit your pull request with a concise explanation of biological and algorithmic rationale.
