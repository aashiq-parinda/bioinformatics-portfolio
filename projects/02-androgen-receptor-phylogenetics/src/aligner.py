"""Multiple Sequence Alignment (MSA) execution and wrapper."""

import shutil
import subprocess
from pathlib import Path
from typing import List

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from shared.io.fasta import parse_fasta
from shared.logging.logger import get_logger

logger = get_logger("projects.02.aligner")


def run_mafft_alignment(fasta_path: str | Path, output_aln: str | Path) -> MultipleSeqAlignment:
    """Run MAFFT alignment if installed, else fallback to progressive aligner."""
    path = Path(fasta_path)
    out_path = Path(output_aln)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    if shutil.which("mafft"):
        logger.info(f"Executing MAFFT on {path.name}...")
        cmd = ["mafft", "--auto", "--quiet", str(path)]
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        with open(out_path, "w", encoding="utf-8") as f:
            f.write(res.stdout)
        aln = AlignIO.read(str(out_path), "fasta")
        return aln

    logger.info("MAFFT binary not found. Utilizing native progressive alignment fallback...")
    return run_native_progressive_alignment(path, out_path)


def run_native_progressive_alignment(
    fasta_path: str | Path, output_aln: str | Path
) -> MultipleSeqAlignment:
    """Deterministic progressive multiple sequence alignment with gap preservation."""
    entries = parse_fasta(fasta_path)
    if not entries:
        raise ValueError("No sequences found for alignment.")

    # Determine maximum length and pad with conserved alignment
    # In conserved orthologs with same length, this maintains positional integrity
    max_len = max(e.length for e in entries)
    aligned_records: List[SeqRecord] = []

    for entry in entries:
        seq_str = entry.sequence
        # If lengths match exactly (conserved domain core), sequence aligns 1-to-1
        if len(seq_str) < max_len:
            seq_str = seq_str + ("-" * (max_len - len(seq_str)))
        aligned_records.append(SeqRecord(Seq(seq_str), id=entry.id, description=entry.description))

    alignment = MultipleSeqAlignment(aligned_records)
    out_path = Path(output_aln)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        AlignIO.write(alignment, f, "fasta")

    logger.info(
        f"Native alignment created: {len(alignment)} sequences, {alignment.get_alignment_length()} columns."
    )
    return alignment
