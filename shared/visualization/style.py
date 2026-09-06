"""Publication-quality matplotlib and seaborn styling setup."""

import matplotlib.pyplot as plt
import seaborn as sns


def apply_scientific_style() -> None:
    """Apply consistent, high-DPI publication styling across all figures."""
    sns.set_theme(style="ticks", font="sans-serif")

    plt.rcParams.update(
        {
            "figure.dpi": 300,
            "savefig.dpi": 300,
            "font.size": 10,
            "axes.labelsize": 11,
            "axes.titlesize": 12,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "figure.titlesize": 14,
            "lines.linewidth": 1.5,
            "axes.linewidth": 1.0,
            "grid.color": "#e0e0e0",
            "grid.linestyle": "--",
            "grid.linewidth": 0.5,
            "figure.autolayout": True,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.1,
        }
    )


# Colorblind-safe palette definitions
PALETTES = {
    "primary": ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"],
    "plddt": {
        "Very High (>90)": "#0053D6",
        "Confident (70-90)": "#65CBF3",
        "Low (50-70)": "#FFDB13",
        "Very Low (<50)": "#FF7D45",
    },
    "diverging": "vlag",
}
