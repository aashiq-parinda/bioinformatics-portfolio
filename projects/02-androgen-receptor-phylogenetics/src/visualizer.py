"""Scientific visualization routines for phylogenetic trees and conservation profiles."""

from pathlib import Path
from typing import Dict, Optional

import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo

from shared.visualization.style import apply_scientific_style


def plot_phylogenetic_tree(
    tree: Phylo.BaseTree.Tree,
    output_path: str | Path,
    title: str = "Vertebrate Androgen Receptor Phylogenetic Tree (NJ)",
) -> Path:
    """Generate a clean, high-DPI phylogram visualization with taxon labels."""
    apply_scientific_style()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 6))

    # Format taxon display names (replace underscores with spaces)
    for clade in tree.find_clades():
        if clade.name:
            clade.name = clade.name.split("|")[0].replace("_", " ")

    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=False,
        label_func=lambda c: c.name if c.is_terminal() else "",
    )

    ax.set_title(title, fontsize=13, fontweight="bold", pad=15)
    ax.set_xlabel("Evolutionary Distance (Substitutions / Site)", fontsize=10)
    ax.set_ylabel("Taxa", fontsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.savefig(out)
    plt.close(fig)
    return out


def plot_conservation_profile(
    conservation_df: pd.DataFrame,
    output_path: str | Path,
    domain_annotation: Optional[Dict[str, Dict]] = None,
    title: str = "Shannon Entropy Conservation Profile Across Alignment",
) -> Path:
    """Plot per-column Shannon entropy with functional domain shading."""
    apply_scientific_style()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(12, 4.5))

    x = conservation_df["position"]
    y = conservation_df["shannon_entropy"]

    # Rolling mean for smooth trendline
    rolling_mean = y.rolling(window=7, min_periods=1, center=True).mean()

    ax.bar(x, y, color="#94a3b8", alpha=0.5, width=1.0, label="Column Entropy (bits)")
    ax.plot(x, rolling_mean, color="#0f766e", lw=1.8, label="Rolling Mean (window=7)")

    # Mark domain boundaries
    if domain_annotation:
        domain_colors = ["#e0f2fe", "#fef3c7", "#dcfce7", "#f3e8ff"]
        for idx, (d_name, d_info) in enumerate(domain_annotation.items()):
            color = domain_colors[idx % len(domain_colors)]
            start = d_info["start_col"]
            end = d_info["end_col"]
            ax.axvspan(start, end, alpha=0.3, color=color, label=f"{d_name} ({start}-{end})")

    ax.set_xlabel("Alignment Position (Residue Column)")
    ax.set_ylabel("Shannon Entropy (bits)")
    ax.set_title(title, fontweight="bold")
    ax.set_ylim(0, max(y.max() * 1.15, 1.0))
    ax.legend(loc="upper right", framealpha=0.9, fontsize=8)

    plt.savefig(out)
    plt.close(fig)
    return out


def plot_domain_conservation(
    domain_entropies: Dict[str, float],
    output_path: str | Path,
    title: str = "Domain-Level Sequence Entropy",
) -> Path:
    """Bar chart comparing average entropy across functional protein domains."""
    apply_scientific_style()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    domains = list(domain_entropies.keys())
    entropies = list(domain_entropies.values())

    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(domains, entropies, color="#1e3a8a", alpha=0.85, edgecolor="#0f172a")

    ax.set_ylabel("Mean Shannon Entropy (bits)")
    ax.set_xlabel("Receptor Domain")
    ax.set_title(title, fontweight="bold")
    ax.set_ylim(0, max(entropies) * 1.25 if entropies else 1.0)

    for bar in bars:
        h = bar.get_height()
        ax.annotate(
            f"{h:.3f}",
            xy=(bar.get_x() + bar.get_width() / 2, h),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=9,
            fontweight="bold",
        )

    plt.savefig(out)
    plt.close(fig)
    return out
