"""Visualizations for AlphaFold structure assessment and confidence profiles."""

from pathlib import Path
from typing import Any, Dict, Optional

import matplotlib.pyplot as plt
import pandas as pd

from shared.visualization.style import apply_scientific_style


def plot_plddt_residue_profile(
    df: pd.DataFrame,
    output_path: str | Path,
    domain_definitions: Optional[Dict[str, Dict[str, Any]]] = None,
    title: str = "AlphaFold Per-Residue Confidence Profile (pLDDT)",
) -> Path:
    """Generate per-residue pLDDT line chart with shaded confidence tiers and domain spans."""
    apply_scientific_style()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(12, 4.5))

    x = df["residue_number"]
    y = df["plddt"]

    # Background confidence bands
    ax.axhspan(90, 100, color="#0053D6", alpha=0.15, label="Very High (>90)")
    ax.axhspan(70, 90, color="#65CBF3", alpha=0.15, label="Confident (70-90)")
    ax.axhspan(50, 70, color="#FFDB13", alpha=0.15, label="Low (50-70)")
    ax.axhspan(0, 50, color="#FF7D45", alpha=0.15, label="Very Low / IDR (<50)")

    # Plot pLDDT line
    ax.plot(x, y, color="#1e293b", lw=1.2)

    # Domain annotations on top
    if domain_definitions:
        domain_colors = ["#cbd5e1", "#bae6fd", "#fef08a", "#bbf7d0"]
        for idx, (d_name, d_info) in enumerate(domain_definitions.items()):
            color = domain_colors[idx % len(domain_colors)]
            start, end = d_info["start"], d_info["end"]
            ax.axvspan(start, end, ymin=0.92, ymax=1.0, color=color, alpha=0.9)
            ax.text(
                (start + end) / 2,
                95,
                d_name,
                ha="center",
                va="center",
                fontsize=8,
                fontweight="bold",
            )

    ax.set_xlabel("Residue Number (UniProt P10275)")
    ax.set_ylabel("Predicted LDDT (0-100)")
    ax.set_title(title, fontweight="bold")
    ax.set_ylim(0, 105)
    ax.set_xlim(x.min() - 5, x.max() + 5)
    ax.legend(loc="lower right", framealpha=0.9, fontsize=8)

    plt.savefig(out)
    plt.close(fig)
    return out


def plot_confidence_distribution(
    summary: Any,
    output_path: str | Path,
    title: str = "AlphaFold Confidence Tier Breakdown",
) -> Path:
    """Bar chart of residue percentages across the four standard confidence tiers."""
    apply_scientific_style()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    tiers = ["Very High (>90)", "Confident (70-90)", "Low (50-70)", "Very Low / IDR (<50)"]
    pcts = [summary.pct_very_high, summary.pct_confident, summary.pct_low, summary.pct_very_low]
    colors = ["#0053D6", "#65CBF3", "#FFDB13", "#FF7D45"]

    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(tiers, pcts, color=colors, edgecolor="#0f172a", alpha=0.85)

    ax.set_ylabel("Residues (%)")
    ax.set_title(title, fontweight="bold")
    ax.set_ylim(0, max(pcts) * 1.2 if pcts else 100)

    for bar in bars:
        h = bar.get_height()
        ax.annotate(
            f"{h:.1f}%",
            xy=(bar.get_x() + bar.get_width() / 2, h),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
        )

    plt.xticks(rotation=15, ha="right")
    plt.savefig(out)
    plt.close(fig)
    return out


def plot_pocket_residues_plddt(
    pocket_df: pd.DataFrame,
    output_path: str | Path,
    title: str = "Active Pocket Residues & AlphaFold Confidence",
) -> Path:
    """Bar chart showing AlphaFold confidence for residues bordering the steroid binding pocket."""
    apply_scientific_style()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    if pocket_df.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No pocket residues available", ha="center", va="center")
        plt.savefig(out)
        plt.close(fig)
        return out

    # Sort by residue number
    sorted_df = pocket_df.sort_values("residue_number").head(15)
    labels = [f"{row['residue_name']}{row['residue_number']}" for _, row in sorted_df.iterrows()]
    scores = sorted_df["plddt"].values

    fig, ax = plt.subplots(figsize=(9, 4))
    colors = ["#0053D6" if s >= 90 else "#65CBF3" if s >= 70 else "#FFDB13" for s in scores]
    ax.bar(labels, scores, color=colors, edgecolor="#0f172a")

    ax.axhline(90, color="#0053D6", linestyle="--", lw=0.8, label="Very High Cutoff (90)")
    ax.set_ylabel("AlphaFold pLDDT Score")
    ax.set_xlabel("Binding Pocket Residue")
    ax.set_title(title, fontweight="bold")
    ax.set_ylim(60, 102)
    ax.legend(loc="lower right")

    plt.xticks(rotation=45, ha="right")
    plt.savefig(out)
    plt.close(fig)
    return out
