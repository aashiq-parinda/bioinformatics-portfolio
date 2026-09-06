"""Unit tests for differential expression and Wald testing."""

import sys
from pathlib import Path

import pandas as pd

project_root = Path(__file__).resolve().parent.parent
repo_root = project_root.parent.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from src.differential_expression import run_differential_expression


def test_run_differential_expression() -> None:
    # 3 control vs 3 treated samples with 1 strongly upregulated gene and 1 invariant
    df = pd.DataFrame(
        {
            "gene_id": ["UP_GENE", "INVARIANT_GENE"],
            "C1": [10.0, 50.0],
            "C2": [12.0, 48.0],
            "C3": [11.0, 52.0],
            "T1": [100.0, 51.0],
            "T2": [105.0, 49.0],
            "T3": [98.0, 50.0],
        }
    )

    res_df, summary = run_differential_expression(
        df,
        control_samples=["C1", "C2", "C3"],
        treated_samples=["T1", "T2", "T3"],
        fdr_cutoff=0.05,
        log2fc_cutoff=1.0,
    )

    assert len(res_df) == 2
    up_row = res_df[res_df["gene_id"] == "UP_GENE"].iloc[0]
    inv_row = res_df[res_df["gene_id"] == "INVARIANT_GENE"].iloc[0]

    assert up_row["status"] == "Up-regulated"
    assert up_row["log2FoldChange"] > 2.0
    assert up_row["padj"] < 0.05

    assert inv_row["status"] == "Not Significant"
    assert abs(inv_row["log2FoldChange"]) < 0.5
