"""Visualizations for molecular docking results and interaction profiling."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from shared.visualization.style import apply_scientific_style


def plot_binding_affinities(
    affinities_df: pd.DataFrame,
    output_path: str | Path,
    title: str = "Molecular Docking Predicted Binding Affinities (ΔG)",
) -> Path:
    """Bar chart comparing binding free energies (kcal/mol) across ligands."""
    apply_scientific_style()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 4.5))

    # Top mode affinity
    top_modes = affinities_df[affinities_df["mode"] == 1].sort_values("affinity_kcal_mol")
    ligands = top_modes["ligand_id"].tolist()
    energies = top_modes["affinity_kcal_mol"].tolist()

    colors = ["#1e3a8a", "#0f766e", "#b45309", "#4338ca"][: len(ligands)]
    bars = ax.bar(ligands, energies, color=colors, edgecolor="#0f172a", alpha=0.85, width=0.55)

    ax.set_ylabel("Binding Affinity ΔG (kcal/mol)")
    ax.set_title(title, fontweight="bold")
    ax.set_ylim(min(energies) * 1.25 if energies else -15, 0)
    ax.axhline(0, color="black", lw=0.8)

    for bar in bars:
        h = bar.get_height()
        ax.annotate(
            f"{h:.2f}",
            xy=(bar.get_x() + bar.get_width() / 2, h),
            xytext=(0, -12),
            textcoords="offset points",
            ha="center",
            va="top",
            color="white",
            fontweight="bold",
            fontsize=10,
        )

    plt.savefig(out)
    plt.close(fig)
    return out


def plot_interaction_heatmap(
    interactions_df: pd.DataFrame,
    output_path: str | Path,
    title: str = "Key Active Site Residue Contact Matrix",
) -> Path:
    """Heatmap showing presence and minimum distance of ligand contacts against active pocket residues."""
    apply_scientific_style()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    if interactions_df.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No interaction contacts detected", ha="center", va="center")
        plt.savefig(out)
        plt.close(fig)
        return out

    # Pivot table: Ligand vs Residue with distance
    pivot = interactions_df.pivot_table(
        index="ligand_id",
        columns="residue_number",
        values="distance_angstroms",
        aggfunc="min",
    )

    fig, ax = plt.subplots(figsize=(10, 4))
    sns.heatmap(
        pivot,
        cmap="mako_r",
        annot=True,
        fmt=".1f",
        cbar_kws={"label": "Minimum Distance (Å)"},
        linewidths=0.5,
        ax=ax,
    )

    ax.set_title(title, fontweight="bold")
    ax.set_xlabel("Receptor Residue Number (Human AR LBD)")
    ax.set_ylabel("Ligand")

    plt.savefig(out)
    plt.close(fig)
    return out


def plot_pose_energy_spectrum(
    all_poses_df: pd.DataFrame,
    output_path: str | Path,
    title: str = "Pose Rank vs. Binding Free Energy Spectrum",
) -> Path:
    """Line plot showing energy spectrum across the 9 conformational modes for each ligand."""
    apply_scientific_style()
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(8, 4))

    for lig_id, grp in all_poses_df.groupby("ligand_id"):
        sorted_grp = grp.sort_values("mode")
        ax.plot(
            sorted_grp["mode"], sorted_grp["affinity_kcal_mol"], marker="o", lw=1.5, label=lig_id
        )

    ax.set_xlabel("Docking Pose Rank (Mode)")
    ax.set_ylabel("Binding Affinity ΔG (kcal/mol)")
    ax.set_title(title, fontweight="bold")
    ax.set_xticks(range(1, 10))
    ax.legend(loc="lower right")

    plt.savefig(out)
    plt.close(fig)
    return out
