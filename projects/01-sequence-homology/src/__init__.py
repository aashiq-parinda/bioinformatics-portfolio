"""Project 01: Automated Protein Sequence & Homology Analysis."""

from .blast import (
    BlastHit,
    filter_blast_hits,
    parse_blast_xml_stream,
    run_local_blast,
    run_remote_ncbi_blast,
)
from .reporter import generate_analysis_report
from .stats import (
    SequenceStatistics,
    calculate_sequence_statistics,
    compute_hydropathy_profile,
)

__all__ = [
    "calculate_sequence_statistics",
    "compute_hydropathy_profile",
    "SequenceStatistics",
    "BlastHit",
    "parse_blast_xml_stream",
    "run_remote_ncbi_blast",
    "run_local_blast",
    "filter_blast_hits",
    "generate_analysis_report",
]
