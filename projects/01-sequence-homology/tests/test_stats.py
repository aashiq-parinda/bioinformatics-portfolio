"""Unit tests for sequence statistics and physicochemical calculations."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from src.stats import (
    calculate_sequence_statistics,
    compute_hydropathy_profile,
)

from shared.utils.validators import validate_protein_sequence


def test_validate_protein_sequence() -> None:
    valid_seq = "MEVQLGLGRVYPRPPSKTYRGAFQNLFQSVREVIQNPGPRHPEAASAAPPGASLLLLQQQ"
    is_valid, msg = validate_protein_sequence(valid_seq)
    assert is_valid
    assert msg == ""

    invalid_seq = "MEVQL123XYZ!"
    is_valid, msg = validate_protein_sequence(invalid_seq, allow_ambiguous=False)
    assert not is_valid
    assert "invalid non-protein characters" in msg


def test_calculate_sequence_statistics_insulin() -> None:
    # Human Insulin sequence (110 aa)
    insulin_seq = "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
    stats = calculate_sequence_statistics(insulin_seq, seq_id="P01308", description="Insulin")

    assert stats.length == 110
    assert 11000 < stats.molecular_weight_da < 13000
    assert 4.5 < stats.isoelectric_point < 7.0
    assert isinstance(stats.gravy_hydropathy, float)
    assert isinstance(stats.amino_acid_composition_pct, dict)
    assert "C" in stats.amino_acid_composition_pct
    assert stats.amino_acid_composition_pct["C"] > 0


def test_compute_hydropathy_profile() -> None:
    test_seq = "MALWMRLLPLLALLALWGPDPAAAFVNQ"
    profile = compute_hydropathy_profile(test_seq, window_size=5)
    assert len(profile) == len(test_seq) - 5 + 1
    assert all(isinstance(score, float) for score in profile)
