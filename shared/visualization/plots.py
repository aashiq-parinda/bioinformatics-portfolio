"""Reusable scientific plotting functions for sequence, structure, and expression figures."""

from pathlib import Path
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .style import apply_scientific_style


def plot_amino_acid_composition(
    composition: Dict[str, float],
    output_path: str | Path,
    title: str = "Amino Acid Composition",
) -> Path:
    """Generate a clean horizontal bar chart of amino acid frequencies (%)."""
    apply_scientific_style()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    sorted_comp = sorted(composition.items(), key=lambda x: x[1], reverse=True)
    aas = [x[0] for x in sorted_comp]
    freqs = [x[1] for x in sorted_comp]

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(aas, freqs, color="#2b5c8f", edgecolor="#1b3b5f", alpha=0.85)

    ax.set_ylabel("Frequency (%)")
    ax.set_xlabel("Amino Acid Residue")
    ax.set_title(title, fontweight="bold")
    ax.set_ylim(0, max(freqs) * 1.15 if freqs else 10)

    for bar in bars:
        height = bar.get_height()
        ax.annotate(
            f"{height:.1f}%",
            xy=(bar.get_x() + bar.get_width() / 2, height),
            xytext=(0, 3),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=7,
        )

    plt.savefig(out)
    plt.close(fig)
    return out


def plot_hydropathy_profile(
    hydropathy_scores: List[float],
    window_size: int,
    output_path: str | Path,
    title: str = "Kyte-Doolittle Hydropathy Profile",
) -> Path:
    """Generate a sliding window Kyte-Doolittle hydropathy profile plot."""
    apply_scientific_style()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 4))
    x = np.arange(1, len(hydropathy_scores) + 1)
    y = np.array(hydropathy_scores)

    ax.plot(x, y, color="#2b5c8f", lw=1.2, label=f"Window Size = {window_size}")
    ax.axhline(0, color="gray", linestyle="--", lw=0.8)
    ax.fill_between(x, y, 0, where=(y >= 0), color="#e27d60", alpha=0.3, label="Hydrophobic")
    ax.fill_between(x, y, 0, where=(y < 0), color="#41b3a3", alpha=0.3, label="Hydrophilic")

    ax.set_xlabel("Residue Position (Center of Window)")
    ax.set_ylabel("Hydropathy Score")
    ax.set_title(title, fontweight="bold")
    ax.legend(loc="upper right", frameon=True)

    plt.savefig(out)
    plt.close(fig)
    return out


def plot_blast_hits_distribution(
    hits_df: pd.DataFrame,
    output_path: str | Path,
    title: str = "Top BLAST Homology Hits",
) -> Path:
    """Generate a scatter/bubble plot of Identity % vs -log10(E-value)."""
    apply_scientific_style()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    if hits_df.empty:
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.text(0.5, 0.5, "No BLAST hits matching filter criteria", ha="center", va="center")
        plt.savefig(out)
        plt.close(fig)
        return out

    fig, ax = plt.subplots(figsize=(8, 5))
    e_values = hits_df["evalue"].clip(lower=1e-180)
    neg_log_eval = -np.log10(e_values.astype(float) + 1e-250)

    scatter = ax.scatter(
        hits_df["identity_pct"],
        neg_log_eval,
        c=hits_df["alignment_length"],
        cmap="viridis",
        s=80,
        alpha=0.85,
        edgecolor="k",
        linewidth=0.5,
    )

    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label("Alignment Length (aa)")

    ax.set_xlabel("Sequence Identity (%)")
    ax.set_ylabel("-log10(E-value)")
    ax.set_title(title, fontweight="bold")

    plt.savefig(out)
    plt.close(fig)
    return out
