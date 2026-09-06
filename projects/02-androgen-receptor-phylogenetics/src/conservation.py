"""Shannon entropy and positional conservation scoring for multiple sequence alignments."""

import math
from collections import Counter
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from Bio.Align import MultipleSeqAlignment

from shared.logging.logger import get_logger

logger = get_logger("projects.02.conservation")


@dataclass
class AlignmentStatistics:
    """Summary metrics of an aligned sequence dataset."""

    num_sequences: int
    alignment_length: int
    gap_percentage: float
    invariant_columns_count: int
    invariant_percentage: float
    mean_entropy: float
    median_entropy: float


def calculate_shannon_entropy(column: str, ignore_gaps: bool = True) -> float:
    """Calculate Shannon entropy for an alignment column in bits.

    H = - sum(p_i * log2(p_i))
    H = 0 implies 100% strict identity/conservation across all taxa.
    """
    chars = [c.upper() for c in column if not (ignore_gaps and c in "-.")]
    if not chars:
        return 0.0

    total = len(chars)
    counts = Counter(chars)
    entropy = 0.0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p)

    return round(entropy, 4)


def analyze_alignment_conservation(
    alignment: MultipleSeqAlignment,
    domain_annotation: Optional[Dict[str, Dict[str, Any]]] = None,
) -> tuple[pd.DataFrame, AlignmentStatistics, Dict[str, float]]:
    """Compute per-position entropy and domain-specific conservation scores.

    Args:
        alignment: Biopython MultipleSeqAlignment object.
        domain_annotation: Optional dictionary defining domains and start/end coordinates.

    Returns:
        Tuple of (per_position_df, alignment_statistics, domain_mean_entropies).
    """
    aln_len = alignment.get_alignment_length()
    num_seqs = len(alignment)

    records: List[Dict[str, Any]] = []
    entropies: List[float] = []
    total_gaps = 0
    total_chars = aln_len * num_seqs

    for col_idx in range(aln_len):
        col_chars = alignment[:, col_idx]
        total_gaps += col_chars.count("-") + col_chars.count(".")
        ent = calculate_shannon_entropy(col_chars, ignore_gaps=True)
        entropies.append(ent)

        # Consensus residue
        counts = Counter([c for c in col_chars if c not in "-."])
        consensus = counts.most_common(1)[0][0] if counts else "-"
        consensus_freq = (
            (counts[consensus] / len([c for c in col_chars if c not in "-."])) if counts else 0.0
        )

        # Determine domain
        assigned_domain = "Unannotated"
        if domain_annotation:
            col_pos = col_idx + 1
            for d_name, d_info in domain_annotation.items():
                if d_info["start_col"] <= col_pos <= d_info["end_col"]:
                    assigned_domain = d_name
                    break

        records.append(
            {
                "position": col_idx + 1,
                "consensus_residue": consensus,
                "consensus_frequency": round(consensus_freq, 3),
                "shannon_entropy": ent,
                "is_invariant": ent == 0.0,
                "domain": assigned_domain,
            }
        )

    df = pd.DataFrame(records)
    inv_count = int(df["is_invariant"].sum())
    stats = AlignmentStatistics(
        num_sequences=num_seqs,
        alignment_length=aln_len,
        gap_percentage=round((total_gaps / total_chars) * 100.0, 2),
        invariant_columns_count=inv_count,
        invariant_percentage=round((inv_count / aln_len) * 100.0, 2),
        mean_entropy=round(float(np.mean(entropies)), 3),
        median_entropy=round(float(np.median(entropies)), 3),
    )

    # Domain summary
    domain_entropies: Dict[str, float] = {}
    for domain, group in df.groupby("domain"):
        domain_entropies[str(domain)] = round(float(group["shannon_entropy"].mean()), 3)

    logger.info(
        f"Conservation analysis: {inv_count}/{aln_len} invariant positions ({stats.invariant_percentage}%), "
        f"Mean Entropy: {stats.mean_entropy:.3f} bits."
    )
    return df, stats, domain_entropies
