"""Scientific visualization library with publication styles and plots."""

from .plots import (
    plot_amino_acid_composition,
    plot_blast_hits_distribution,
    plot_hydropathy_profile,
)
from .style import PALETTES, apply_scientific_style

__all__ = [
    "apply_scientific_style",
    "PALETTES",
    "plot_amino_acid_composition",
    "plot_hydropathy_profile",
    "plot_blast_hits_distribution",
]
