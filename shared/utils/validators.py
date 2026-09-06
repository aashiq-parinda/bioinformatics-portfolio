"""Validation functions for sequences, accession IDs, and file formats."""

import re
from pathlib import Path
from typing import Set

# Standard IUPAC amino acid codes (including Selenocysteine U and Pyrrolysine O, and ambiguity codes)
IUPAC_AMINO_ACIDS: Set[str] = set("ACDEFGHIKLMNPQRSTVWYUO")
IUPAC_AMBIGUOUS_AMINO_ACIDS: Set[str] = set("ACDEFGHIKLMNPQRSTVWYUOBZJX*")


def validate_protein_sequence(sequence: str, allow_ambiguous: bool = True) -> tuple[bool, str]:
    """Validate whether a given string is a valid amino acid sequence.

    Args:
        sequence: The raw amino acid sequence.
        allow_ambiguous: If True, allow B, Z, J, X, * characters.

    Returns:
        Tuple of (is_valid: bool, error_message: str).
    """
    clean_seq = "".join(sequence.split()).upper()
    if not clean_seq:
        return False, "Sequence is empty."

    valid_alphabet = IUPAC_AMBIGUOUS_AMINO_ACIDS if allow_ambiguous else IUPAC_AMINO_ACIDS
    invalid_chars = set(clean_seq) - valid_alphabet

    if invalid_chars:
        return (
            False,
            f"Sequence contains invalid non-protein characters: {', '.join(sorted(invalid_chars))}",
        )

    return True, ""


def validate_file_exists(file_path: str | Path, must_not_be_empty: bool = True) -> tuple[bool, str]:
    """Verify that a given file exists and optionally is non-empty."""
    path = Path(file_path)
    if not path.is_file():
        return False, f"File does not exist: {file_path}"
    if must_not_be_empty and path.stat().st_size == 0:
        return False, f"File is empty: {file_path}"
    return True, ""


def sanitize_filename(name: str) -> str:
    """Sanitize strings for safe file system usage."""
    return re.sub(r"[^a-zA-Z0-9_\-\.]", "_", name)
