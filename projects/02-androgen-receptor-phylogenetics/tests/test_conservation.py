"""Unit tests for Shannon entropy and alignment conservation analysis."""

import math
import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
repo_root = project_root.parent.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from src.conservation import analyze_alignment_conservation, calculate_shannon_entropy


def test_calculate_shannon_entropy_identical() -> None:
    # 100% identical column must yield H = 0.0
    col = "AAAAAAAAAA"
    h = calculate_shannon_entropy(col)
    assert h == 0.0


def test_calculate_shannon_entropy_two_state() -> None:
    # 50/50 two-state distribution must yield H = 1.0 bit
    col = "AAAAABBBBB"
    h = calculate_shannon_entropy(col)
    assert math.isclose(h, 1.0, abs_tol=1e-3)


def test_analyze_alignment_conservation() -> None:
    seqs = [
        SeqRecord(Seq("ACDEF"), id="Seq1"),
        SeqRecord(Seq("ACDEG"), id="Seq2"),
        SeqRecord(Seq("ACDEF"), id="Seq3"),
    ]
    alignment = MultipleSeqAlignment(seqs)
    domain_annot = {
        "DomainA": {"start_col": 1, "end_col": 3},
        "DomainB": {"start_col": 4, "end_col": 5},
    }

    df, stats, dom_entropies = analyze_alignment_conservation(alignment, domain_annot)

    assert stats.num_sequences == 3
    assert stats.alignment_length == 5
    assert stats.invariant_columns_count == 4  # A, C, D, E are invariant
    assert "DomainA" in dom_entropies
    assert dom_entropies["DomainA"] == 0.0  # All positions 1-3 are invariant
