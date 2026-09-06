"""Unit tests for DESeq2 Median-of-Ratios count normalization."""

import sys
from pathlib import Path

import pandas as pd

project_root = Path(__file__).resolve().parent.parent
repo_root = project_root.parent.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from src.normalization import calculate_size_factors, normalize_counts


def test_calculate_size_factors() -> None:
    df = pd.DataFrame(
        {
            "gene_id": ["G1", "G2", "G3", "G4"],
            "S1": [100, 200, 300, 400],
            "S2": [200, 400, 600, 800],  # Exactly 2x sequencing depth
        }
    )
    samples = ["S1", "S2"]
    factors = calculate_size_factors(df, samples)

    assert "S1" in factors and "S2" in factors
    # S2 size factor should be approximately 2x S1
    assert factors["S2"] > factors["S1"]

    norm_df, log2_df = normalize_counts(df, samples, factors)
    assert "S1" in norm_df.columns
    assert "S2" in norm_df.columns
    assert "S1" in log2_df.columns
