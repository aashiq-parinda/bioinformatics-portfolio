"""DESeq2 Median-of-Ratios count normalization and log2 transformation."""

from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from shared.logging.logger import get_logger

logger = get_logger("projects.05.normalization")


def calculate_size_factors(counts_df: pd.DataFrame, sample_cols: List[str]) -> Dict[str, float]:
    """Compute DESeq2 median-of-ratios size factors for sequencing depth normalization.

    Args:
        counts_df: DataFrame with genes as rows and samples as columns.
        sample_cols: List of sample column names.

    Returns:
        Dictionary mapping sample_id to its computed size factor.
    """
    sub_counts = counts_df[sample_cols].astype(float).copy()

    # Step 1: Geometric mean per gene across all samples (ignoring rows with 0s)
    # Using log-space mean to prevent floating point overflow/underflow
    log_counts = np.log(sub_counts.replace(0, np.nan))
    pseudo_refs = np.exp(log_counts.mean(axis=1))

    # Step 2: Ratio of each sample count to the pseudo-reference
    size_factors: Dict[str, float] = {}
    for sample in sample_cols:
        ratios = sub_counts[sample] / pseudo_refs
        # Step 3: Median of finite non-zero ratios
        valid_ratios = ratios[np.isfinite(ratios) & (ratios > 0)]
        sf = float(np.median(valid_ratios)) if not valid_ratios.empty else 1.0
        size_factors[sample] = round(sf, 4)

    logger.info(f"Computed size factors: {size_factors}")
    return size_factors


def normalize_counts(
    counts_df: pd.DataFrame,
    sample_cols: List[str],
    size_factors: Dict[str, float],
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Normalize raw counts by dividing by size factors, and compute log2-normalized counts.

    Returns:
        Tuple of (normalized_counts_df, log2_normalized_df).
    """
    norm_df = counts_df.copy()
    log2_df = counts_df.copy()

    for sample in sample_cols:
        sf = size_factors.get(sample, 1.0)
        norm_df[sample] = round(norm_df[sample] / sf, 2)
        log2_df[sample] = round(np.log2(norm_df[sample] + 1.0), 3)

    return norm_df, log2_df
