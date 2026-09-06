"""Negative Binomial Wald hypothesis testing and differential expression statistics."""

from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from shared.logging.logger import get_logger

logger = get_logger("projects.05.de")


@dataclass
class DESummary:
    """High-level summary of differential expression testing."""

    total_genes_tested: int
    num_significant_up: int
    num_significant_down: int
    num_not_significant: int
    fdr_threshold: float
    log2fc_threshold: float


def benjamini_hochberg_fdr(p_values: np.ndarray) -> np.ndarray:
    """Compute Benjamini-Hochberg False Discovery Rate (FDR) adjusted p-values."""
    n = len(p_values)
    if n == 0:
        return np.array([])

    sorted_indices = np.argsort(p_values)
    sorted_p = p_values[sorted_indices]
    fdr = np.empty(n, dtype=float)

    # Calculate cumulative minimum of (p * n / rank)
    cum_min = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        val = (sorted_p[i] * n) / rank
        cum_min = min(cum_min, val)
        fdr[i] = min(cum_min, 1.0)

    # Reorder to original indices
    adj_p = np.empty(n, dtype=float)
    adj_p[sorted_indices] = fdr
    return adj_p


def run_differential_expression(
    norm_counts_df: pd.DataFrame,
    control_samples: List[str],
    treated_samples: List[str],
    fdr_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
    min_count_sum: int = 10,
) -> Tuple[pd.DataFrame, DESummary]:
    """Execute Negative Binomial Wald test for two-group comparative transcriptomics.

    Args:
        norm_counts_df: DataFrame of normalized counts with 'gene_id' column.
        control_samples: Column names for reference group.
        treated_samples: Column names for treatment group.
        fdr_cutoff: Significance threshold for FDR-adjusted p-value.
        log2fc_cutoff: Absolute log2 fold change threshold.
        min_count_sum: Minimum sum of counts across all samples to retain gene.

    Returns:
        Tuple of (results_dataframe, DESummary).
    """
    all_samples = control_samples + treated_samples
    df = norm_counts_df.copy()

    # Filter lowly expressed genes
    total_counts = df[all_samples].sum(axis=1)
    df = df[total_counts >= min_count_sum].reset_index(drop=True)

    n_ctrl = len(control_samples)
    n_treat = len(treated_samples)

    ctrl_means = df[control_samples].mean(axis=1).values
    treat_means = df[treated_samples].mean(axis=1).values
    base_means = df[all_samples].mean(axis=1).values

    # Compute Log2 Fold Change with pseudocount
    log2_fc = np.log2((treat_means + 0.5) / (ctrl_means + 0.5))

    # Dispersion estimation (empirical Bayes approximation)
    ctrl_var = df[control_samples].var(axis=1).values
    treat_var = df[treated_samples].var(axis=1).values
    pooled_var = (ctrl_var + treat_var) / 2.0
    # Overdispersion alpha = max(0, (var - mu) / mu^2)
    alpha = np.clip((pooled_var - base_means) / (base_means**2 + 1e-6), 0.01, 2.0)

    # Standard error of log2 fold change
    se_ctrl = 1.0 / (np.sum(df[control_samples].values, axis=1) + 1e-6) + alpha / n_ctrl
    se_treat = 1.0 / (np.sum(df[treated_samples].values, axis=1) + 1e-6) + alpha / n_treat
    lfc_se = np.sqrt(np.clip(se_ctrl + se_treat, 1e-6, 10.0)) / np.log(2.0)

    # Wald statistic and p-value
    wald_stat = log2_fc / (lfc_se + 1e-6)
    p_values = 2.0 * stats.norm.sf(np.abs(wald_stat))
    p_values = np.clip(p_values, 1e-300, 1.0)

    # FDR adjustment
    adj_p = benjamini_hochberg_fdr(p_values)

    # Assign Significance Status
    status = []
    for lfc, fdr in zip(log2_fc, adj_p):
        if fdr <= fdr_cutoff:
            if lfc >= log2fc_cutoff:
                status.append("Up-regulated")
            elif lfc <= -log2fc_cutoff:
                status.append("Down-regulated")
            else:
                status.append("Not Significant")
        else:
            status.append("Not Significant")

    res_df = pd.DataFrame(
        {
            "gene_id": df["gene_id"],
            "baseMean": np.round(base_means, 2),
            "log2FoldChange": np.round(log2_fc, 3),
            "lfcSE": np.round(lfc_se, 3),
            "stat": np.round(wald_stat, 2),
            "pvalue": p_values,
            "padj": adj_p,
            "status": status,
        }
    )

    # Sort by padj ascending, then abs(log2FC) descending
    res_df = res_df.sort_values(["padj", "pvalue"]).reset_index(drop=True)

    n_up = int(np.sum(res_df["status"] == "Up-regulated"))
    n_down = int(np.sum(res_df["status"] == "Down-regulated"))
    n_ns = int(np.sum(res_df["status"] == "Not Significant"))

    summary = DESummary(
        total_genes_tested=len(res_df),
        num_significant_up=n_up,
        num_significant_down=n_down,
        num_not_significant=n_ns,
        fdr_threshold=fdr_cutoff,
        log2fc_threshold=log2fc_cutoff,
    )

    logger.info(
        f"Differential Expression: {len(res_df)} genes tested. "
        f"Significant Up: {n_up}, Significant Down: {n_down}, Non-significant: {n_ns}."
    )
    return res_df, summary
