"""BLAST execution, hit parsing, and homology filtering engine."""

import os
import shutil
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd
from Bio.Blast import NCBIWWW, NCBIXML

from shared.logging.logger import get_logger
from shared.utils.rate_limiter import RateLimiter, retry_with_backoff

logger = get_logger("projects.01.blast")
ncbi_limiter = RateLimiter(max_calls_per_second=2.0)


@dataclass
class BlastHit:
    """Standardized representation of a single BLAST homology hit."""

    query_id: str
    hit_id: str
    accession: str
    title: str
    organism: str
    evalue: float
    bit_score: float
    raw_score: float
    identity_pct: float
    positive_pct: float
    query_coverage_pct: float
    alignment_length: int
    query_start: int
    query_end: int
    hit_start: int
    hit_end: int

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def parse_blast_xml_stream(xml_handle: Any, query_length: int) -> List[BlastHit]:
    """Parse BLAST XML output stream and extract standardized BlastHit records."""
    hits: List[BlastHit] = []
    try:
        blast_records = NCBIXML.parse(xml_handle)
        for record in blast_records:
            q_len = query_length or record.query_length or 1
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    ident_pct = round((hsp.identities / hsp.align_length) * 100.0, 2)
                    pos_pct = round((hsp.positives / hsp.align_length) * 100.0, 2)
                    q_cov = round((abs(hsp.query_end - hsp.query_start + 1) / q_len) * 100.0, 2)

                    # Extract accession and organism from title
                    title_str = alignment.title
                    accession_str = alignment.accession or alignment.hit_id
                    organism_str = "Unknown"
                    if "[" in title_str and "]" in title_str:
                        organism_str = title_str.split("[")[-1].split("]")[0]

                    hits.append(
                        BlastHit(
                            query_id=record.query_id or "Query",
                            hit_id=alignment.hit_id,
                            accession=accession_str,
                            title=title_str,
                            organism=organism_str,
                            evalue=float(hsp.expect),
                            bit_score=float(hsp.bits),
                            raw_score=float(hsp.score),
                            identity_pct=ident_pct,
                            positive_pct=pos_pct,
                            query_coverage_pct=q_cov,
                            alignment_length=int(hsp.align_length),
                            query_start=int(hsp.query_start),
                            query_end=int(hsp.query_end),
                            hit_start=int(hsp.sbjct_start),
                            hit_end=int(hsp.sbjct_end),
                        )
                    )
    except Exception as e:
        logger.error(f"Error parsing BLAST XML: {e}")
        raise

    return hits


@retry_with_backoff(max_retries=3, base_delay=3.0)
def run_remote_ncbi_blast(
    sequence: str,
    database: str = "swissprot",
    evalue_threshold: float = 0.001,
    hitlist_size: int = 50,
) -> str:
    """Execute remote NCBI QBLAST search with polite rate-limiting."""
    ncbi_limiter.wait()
    logger.info(
        f"Submitting QBLAST query to NCBI (database={database}, evalue={evalue_threshold})..."
    )

    result_handle = NCBIWWW.qblast(
        program="blastp",
        database=database,
        sequence=sequence,
        expect=evalue_threshold,
        hitlist_size=hitlist_size,
    )
    xml_data = result_handle.read()
    result_handle.close()
    logger.info(f"Received QBLAST response ({len(xml_data)} bytes)")
    return xml_data


def run_local_blast(
    query_fasta: str | Path,
    database_path: str,
    evalue: float = 0.001,
    max_target_seqs: int = 50,
    num_threads: int = 4,
) -> str:
    """Execute local blastp command line if installed."""
    if not shutil.which("blastp"):
        raise FileNotFoundError("Local blastp executable not found in PATH.")

    with tempfile.NamedTemporaryFile(suffix=".xml", delete=False) as tmp_out:
        tmp_out_path = tmp_out.name

    cmd = [
        "blastp",
        "-query",
        str(query_fasta),
        "-db",
        database_path,
        "-evalue",
        str(evalue),
        "-max_target_seqs",
        str(max_target_seqs),
        "-num_threads",
        str(num_threads),
        "-outfmt",
        "5",  # XML format
        "-out",
        tmp_out_path,
    ]

    logger.info(f"Running local BLAST: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

    with open(tmp_out_path, "r", encoding="utf-8") as f:
        xml_content = f.read()

    os.remove(tmp_out_path)
    return xml_content


def filter_blast_hits(
    hits: List[BlastHit],
    min_identity: float = 30.0,
    max_evalue: float = 1e-5,
    min_query_coverage: float = 40.0,
    max_hits: int = 20,
) -> List[BlastHit]:
    """Filter BLAST hits based on identity, E-value, and coverage thresholds."""
    filtered: List[BlastHit] = []
    for hit in hits:
        if hit.evalue > max_evalue:
            continue
        if hit.identity_pct < min_identity:
            continue
        if hit.query_coverage_pct < min_query_coverage:
            continue
        filtered.append(hit)

    # Sort primarily by E-value ascending, then bit score descending
    filtered.sort(key=lambda h: (h.evalue, -h.bit_score))
    return filtered[:max_hits]


def hits_to_dataframe(hits: List[BlastHit]) -> pd.DataFrame:
    """Convert a list of BlastHit objects into a pandas DataFrame."""
    if not hits:
        return pd.DataFrame(
            columns=[
                "accession",
                "organism",
                "identity_pct",
                "positive_pct",
                "query_coverage_pct",
                "evalue",
                "bit_score",
                "alignment_length",
            ]
        )
    return pd.DataFrame([h.to_dict() for h in hits])
