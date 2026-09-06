"""Robust FASTA reader, writer, and validator using Biopython."""

from dataclasses import dataclass
from pathlib import Path
from typing import List

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from shared.logging.logger import get_logger
from shared.utils.validators import validate_protein_sequence

logger = get_logger("shared.io.fasta")


@dataclass
class SequenceEntry:
    """Structured representation of a biological sequence."""

    id: str
    description: str
    sequence: str

    @property
    def length(self) -> int:
        return len(self.sequence)

    def to_seq_record(self) -> SeqRecord:
        return SeqRecord(Seq(self.sequence), id=self.id, description=self.description)


def parse_fasta(file_path: str | Path, validate_protein: bool = True) -> List[SequenceEntry]:
    """Parse a FASTA file and return a list of SequenceEntry objects.

    Args:
        file_path: Path to the FASTA file.
        validate_protein: If True, validates against protein alphabet.

    Returns:
        List of SequenceEntry dataclass instances.
    """
    path = Path(file_path)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    entries: List[SequenceEntry] = []
    for record in SeqIO.parse(str(path), "fasta"):
        seq_str = str(record.seq).strip().upper()
        if validate_protein:
            is_valid, msg = validate_protein_sequence(seq_str)
            if not is_valid:
                logger.warning(f"Record {record.id} in {path.name}: {msg}")
        entries.append(
            SequenceEntry(
                id=record.id,
                description=record.description,
                sequence=seq_str,
            )
        )

    if not entries:
        raise ValueError(f"No valid FASTA records found in {path}")

    logger.info(f"Loaded {len(entries)} sequences from {path.name}")
    return entries


def write_fasta(entries: List[SequenceEntry] | List[SeqRecord], output_path: str | Path) -> None:
    """Write sequence records to a FASTA file."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    records: List[SeqRecord] = []
    for entry in entries:
        if isinstance(entry, SequenceEntry):
            records.append(entry.to_seq_record())
        else:
            records.append(entry)

    with open(path, "w", encoding="utf-8") as f:
        SeqIO.write(records, f, "fasta")
    logger.info(f"Wrote {len(records)} records to {path}")
