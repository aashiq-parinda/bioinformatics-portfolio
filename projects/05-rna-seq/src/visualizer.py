"""Publication-grade visualizations for transcriptomics, PCA, Volcano, and Heatmaps."""

from pathlib import Path
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from shared.visualization.style import apply_scientific_style


def plot_pca(
    log2_counts_df: pd.DataFrame,
    sample_cols: List[str],
    metadata_df: pd.DataFrame,
    output_path: str | Path,
    title: str = "Principal Component Analysis (PCA)",
) -> Path:
    """Compute and plot PCA of samples based on top 500 most variable genes."""
    apply_scientific_style()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    matrix = log2_counts_df[sample_cols].values
    variances = np.var(matrix, axis=1)
    top_indices = np.argsort(variances)[::-1][:500]
    top_matrix = matrix[top_indices, :].T  # shape: (n_samples, n_genes)

    # Standardize / center matrix across samples
    centered = top_matrix - np.mean(top_matrix, axis=0)
    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    pcs = u * s
    var_exp = (s**2) / np.sum(s**2) * 100.0

    meta_map = dict(zip(metadata_df["sample_id"], metadata_df["condition"]))
    conditions = [meta_map.get(s, "Unknown") for s in sample_cols]

    fig, ax = plt.subplots(figsize=(7, 5))
    color_map = {"Control": "#2563eb", "Treated": "#dc2626"}

    for s_idx, sample in enumerate(sample_cols):
        cond = conditions[s_idx]
        color = color_map.get(cond, "#64748b")
        ax.scatter(
            pcs[s_idx, 0],
            pcs[s_idx, 1],
            c=color,
            s=120,
            edgecolors="k",
            linewidth=0.8,
            label=cond if cond not in ax.get_legend_handles_labels()[1] else "",
        )
        ax.annotate(
            sample,
            (pcs[s_idx, 0], pcs[s_idx, 1]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )

    ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}% Variance)")
    ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}% Variance)")
    ax.set_title(title, fontweight="bold")
    ax.legend(loc="upper right", frameon=True)

    plt.savefig(out)
    plt.close(fig)
    return out


def plot_volcano(
    de_results_df: pd.DataFrame,
    output_path: str | Path,
    fdr_cutoff: float = 0.05,
    log2fc_cutoff: float = 1.0,
    top_genes_to_label: Optional[List[str]] = None,
    title: str = "Volcano Plot: Differential Expression",
) -> Path:
    """Generate Volcano plot (-log10 FDR vs log2FC) with colored status and labeled top genes."""
    apply_scientific_style()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 6))

    df = de_results_df.copy()
    df["neg_log10_padj"] = -np.log10(df["padj"].clip(lower=1e-100))

    # Scatter categories
    up = df[df["status"] == "Up-regulated"]
    down = df[df["status"] == "Down-regulated"]
    ns = df[df["status"] == "Not Significant"]

    ax.scatter(
        ns["log2FoldChange"],
        ns["neg_log10_padj"],
        c="#94a3b8",
        s=25,
        alpha=0.5,
        label=f"Not Sig ({len(ns)})",
    )
    ax.scatter(
        up["log2FoldChange"],
        up["neg_log10_padj"],
        c="#dc2626",
        s=35,
        alpha=0.85,
        label=f"Up ({len(up)})",
    )
    ax.scatter(
        down["log2FoldChange"],
        down["neg_log10_padj"],
        c="#2563eb",
        s=35,
        alpha=0.85,
        label=f"Down ({len(down)})",
    )

    # Cutoff threshold lines
    ax.axvline(log2fc_cutoff, color="gray", linestyle="--", lw=0.8)
    ax.axvline(-log2fc_cutoff, color="gray", linestyle="--", lw=0.8)
    ax.axhline(-np.log10(fdr_cutoff), color="gray", linestyle="--", lw=0.8)

    # Label top canonical genes
    labels = top_genes_to_label or ["KLK3", "FKBP5", "TMPRSS2", "NKX3-1", "MYC", "SLC45A3"]
    for _, row in df.iterrows():
        if row["gene_id"] in labels:
            ax.annotate(
                row["gene_id"],
                (row["log2FoldChange"], row["neg_log10_padj"]),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=8,
                fontweight="bold",
            )

    ax.set_xlabel("log2(Fold Change)")
    ax.set_ylabel("-log10(FDR Adjusted p-value)")
    ax.set_title(title, fontweight="bold")
    ax.legend(loc="upper left")

    plt.savefig(out)
    plt.close(fig)
    return out


def plot_ma(
    de_results_df: pd.DataFrame,
    output_path: str | Path,
    title: str = "MA Plot: log2FC vs Mean Expression",
) -> Path:
    """Generate MA plot (log2FC vs log10 baseMean)."""
    apply_scientific_style()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 5))

    df = de_results_df.copy()
    sig = df[df["status"] != "Not Significant"]
    ns = df[df["status"] == "Not Significant"]

    ax.scatter(
        np.log10(ns["baseMean"].clip(lower=1.0)),
        ns["log2FoldChange"],
        c="#94a3b8",
        s=20,
        alpha=0.4,
        label="Not Significant",
    )
    ax.scatter(
        np.log10(sig["baseMean"].clip(lower=1.0)),
        sig["log2FoldChange"],
        c="#dc2626",
        s=30,
        alpha=0.85,
        label="Significant DE",
    )

    ax.axhline(0, color="blue", linestyle="-", lw=0.8)
    ax.set_xlabel("log10(Mean Normalized Counts)")
    ax.set_ylabel("log2(Fold Change)")
    ax.set_title(title, fontweight="bold")
    ax.legend(loc="upper right")

    plt.savefig(out)
    plt.close(fig)
    return out


def plot_expression_heatmap(
    log2_counts_df: pd.DataFrame,
    de_results_df: pd.DataFrame,
    sample_cols: List[str],
    output_path: str | Path,
    n_top: int = 25,
    title: str = "Top Differentially Expressed Genes (Z-score)",
) -> Path:
    """Clustered heatmap of top DE genes standardized by z-score."""
    apply_scientific_style()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Pick top n genes by padj
    top_genes = de_results_df.head(n_top)["gene_id"].tolist()
    sub_df = log2_counts_df[log2_counts_df["gene_id"].isin(top_genes)].set_index("gene_id")[
        sample_cols
    ]

    # Row-wise Z-score standardization
    z_scores = sub_df.apply(lambda x: (x - x.mean()) / (x.std() + 1e-6), axis=1)

    fig, ax = plt.subplots(figsize=(8, 8))
    sns.heatmap(
        z_scores,
        cmap="vlag",
        center=0,
        linewidths=0.5,
        cbar_kws={"label": "Z-score"},
        ax=ax,
    )

    ax.set_title(title, fontweight="bold", pad=15)
    ax.set_xlabel("Sequencing Sample")
    ax.set_ylabel("Gene")

    plt.savefig(out)
    plt.close(fig)
    return out


def plot_pathway_enrichment(
    enrichment_df: pd.DataFrame,
    output_path: str | Path,
    title: str = "Enriched Functional Pathways & Gene Ontology",
) -> Path:
    """Horizontal bar chart of top enriched pathways ranked by -log10 p-value."""
    apply_scientific_style()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    if enrichment_df.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No significantly enriched pathways", ha="center", va="center")
        plt.savefig(out)
        plt.close(fig)
        return out

    top_df = enrichment_df.head(10).sort_values("p_value", ascending=True)
    y_labels = [row["pathway_name"][:35] for _, row in top_df.iterrows()]
    x_vals = top_df["neg_log10_p"].values

    fig, ax = plt.subplots(figsize=(9, 4.5))
    ax.barh(y_labels, x_vals, color="#0f766e", alpha=0.85, edgecolor="#042f2e")

    ax.set_xlabel("-log10(p-value)")
    ax.set_title(title, fontweight="bold")
    ax.axvline(-np.log10(0.05), color="red", linestyle="--", lw=0.8, label="p = 0.05 Cutoff")
    ax.legend(loc="lower right")

    plt.savefig(out)
    plt.close(fig)
    return out
