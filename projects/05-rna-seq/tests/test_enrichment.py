"""Unit tests for pathway over-representation analysis."""

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
repo_root = project_root.parent.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from src.pathway_enrichment import run_pathway_enrichment


def test_run_pathway_enrichment() -> None:
    # Significant genes enriched in steroid hormone response
    sig_genes = {"KLK3", "FKBP5", "TMPRSS2", "NKX3-1", "SLC45A3", "PMEPA1"}
    background_genes = sig_genes | {f"GENE_{i}" for i in range(500)}

    res_df = run_pathway_enrichment(sig_genes, background_genes, p_val_cutoff=0.05)

    assert not res_df.empty
    assert "pathway_name" in res_df.columns
    assert "p_value" in res_df.columns
    assert res_df.iloc[0]["p_value"] < 0.05
