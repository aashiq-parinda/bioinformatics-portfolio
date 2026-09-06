"""Physicochemical sequence calculations and amino-acid statistics."""

from dataclasses import asdict, dataclass
from typing import Dict, List

import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# Standard Kyte-Doolittle Hydropathy Scale
KYTE_DOOLITTLE_SCALE: Dict[str, float] = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}


@dataclass
class SequenceStatistics:
    """Computed physicochemical metrics for an amino acid sequence."""

    sequence_id: str
    description: str
    length: int
    molecular_weight_da: float
    isoelectric_point: float
    gravy_hydropathy: float
    aromaticity: float
    instability_index: float
    is_stable: bool
    amino_acid_composition_pct: Dict[str, float]
    extinction_coefficient_reduced: int
    extinction_coefficient_cystines: int

    def to_dict(self) -> Dict:
        return asdict(self)


def calculate_sequence_statistics(
    sequence: str, seq_id: str = "seq", description: str = ""
) -> SequenceStatistics:
    """Calculate comprehensive physicochemical statistics for a protein sequence.

    Args:
        sequence: Clean amino acid sequence string.
        seq_id: Sequence identifier.
        description: Sequence header description.

    Returns:
        SequenceStatistics dataclass.
    """
    clean_seq = "".join(sequence.split()).upper()
    # Filter non-standard characters for strict ProtParam calculations
    standard_aa_only = "".join([aa for aa in clean_seq if aa in KYTE_DOOLITTLE_SCALE])
    if not standard_aa_only:
        raise ValueError(f"Sequence {seq_id} contains no standard amino acids.")

    analysis = ProteinAnalysis(standard_aa_only)

    # Core physicochemical properties
    length = len(clean_seq)
    mw = round(analysis.molecular_weight(), 2)
    pi = round(analysis.isoelectric_point(), 2)
    gravy = round(analysis.gravy(), 3)
    aromaticity = round(analysis.aromaticity(), 3)

    try:
        instability = round(analysis.instability_index(), 2)
        is_stable = instability < 40.0
    except Exception:
        instability = 0.0
        is_stable = True

    # Amino acid percentage composition
    if hasattr(analysis, "amino_acids_percent"):
        raw_comp = analysis.amino_acids_percent
    elif hasattr(analysis, "get_amino_acids_percent"):
        raw_comp = analysis.get_amino_acids_percent()
    else:
        from collections import Counter

        counts = Counter(standard_aa_only)
        raw_comp = {aa: counts.get(aa, 0) / len(standard_aa_only) for aa in KYTE_DOOLITTLE_SCALE}
    comp_pct = {aa: round(pct * 100, 2) for aa, pct in sorted(raw_comp.items())}

    # Extinction coefficients (M^-1 cm^-1 at 280 nm)
    try:
        ext_coeff = analysis.molar_extinction_coefficient()
        ext_red, ext_cys = ext_coeff[0], ext_coeff[1]
    except Exception:
        ext_red, ext_cys = 0, 0

    return SequenceStatistics(
        sequence_id=seq_id,
        description=description,
        length=length,
        molecular_weight_da=mw,
        isoelectric_point=pi,
        gravy_hydropathy=gravy,
        aromaticity=aromaticity,
        instability_index=instability,
        is_stable=is_stable,
        amino_acid_composition_pct=comp_pct,
        extinction_coefficient_reduced=ext_red,
        extinction_coefficient_cystines=ext_cys,
    )


def compute_hydropathy_profile(sequence: str, window_size: int = 9) -> List[float]:
    """Compute a sliding-window Kyte-Doolittle hydropathy profile."""
    clean_seq = "".join(sequence.split()).upper()
    if len(clean_seq) < window_size:
        return [float(np.mean([KYTE_DOOLITTLE_SCALE.get(aa, 0.0) for aa in clean_seq]))]

    scores: List[float] = []

    for i in range(len(clean_seq) - window_size + 1):
        window = clean_seq[i : i + window_size]
        avg_score = sum(KYTE_DOOLITTLE_SCALE.get(aa, 0.0) for aa in window) / window_size
        scores.append(round(avg_score, 3))

    return scores
