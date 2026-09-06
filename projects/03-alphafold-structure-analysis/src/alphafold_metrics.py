"""Extraction and statistical evaluation of AlphaFold pLDDT confidence metrics."""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, is_aa

from shared.logging.logger import get_logger

logger = get_logger("projects.03.plddt")


@dataclass
class ConfidenceSummary:
    """Statistical summary of AlphaFold pLDDT scores."""

    total_residues: int
    mean_plddt: float
    median_plddt: float
    pct_very_high: float  # >90
    pct_confident: float  # 70-90
    pct_low: float  # 50-70
    pct_very_low: float  # <50 (IDR)
    domain_metrics: Dict[str, Dict[str, float]]


def extract_alphafold_plddt(
    pdb_path: str | Path,
    domain_definitions: Optional[Dict[str, Dict[str, Any]]] = None,
) -> tuple[pd.DataFrame, ConfidenceSummary]:
    """Parse an AlphaFold PDB file and extract per-residue pLDDT scores and tier categories.

    Args:
        pdb_path: Path to AlphaFold PDB file.
        domain_definitions: Optional dictionary defining domains and residue spans.

    Returns:
        Tuple of (residue_dataframe, ConfidenceSummary).
    """
    path = Path(pdb_path)
    if not path.is_file():
        raise FileNotFoundError(f"AlphaFold PDB file not found: {path}")

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("alphafold", str(path))

    records: List[Dict[str, Any]] = []
    plddt_list: List[float] = []

    for model in structure:
        for chain in model:
            for res in chain:
                if not is_aa(res, standard=True):
                    continue
                res_id = res.get_id()[1]
                res_name = res.get_resname()

                # AlphaFold stores per-residue pLDDT in B-factor
                ca_atoms = [atom for atom in res if atom.get_name() == "CA"]
                if ca_atoms:
                    score = float(ca_atoms[0].get_bfactor())
                else:
                    b_factors = [atom.get_bfactor() for atom in res.get_atoms()]
                    score = float(np.mean(b_factors)) if b_factors else 0.0

                plddt_list.append(score)

                # Confidence category
                if score >= 90.0:
                    tier = "Very High (>90)"
                elif score >= 70.0:
                    tier = "Confident (70-90)"
                elif score >= 50.0:
                    tier = "Low (50-70)"
                else:
                    tier = "Very Low (<50, IDR)"

                # Assign domain
                domain_name = "Unannotated"
                if domain_definitions:
                    for d_key, d_info in domain_definitions.items():
                        if d_info["start"] <= res_id <= d_info["end"]:
                            domain_name = d_key
                            break

                records.append(
                    {
                        "residue_number": res_id,
                        "residue_name": res_name,
                        "chain": chain.get_id(),
                        "plddt": round(score, 2),
                        "confidence_tier": tier,
                        "domain": domain_name,
                    }
                )

    df = pd.DataFrame(records)
    total = len(df)
    if total == 0:
        raise ValueError("No standard amino acid residues found in PDB.")

    scores = np.array(plddt_list)
    vh = float(np.sum(scores >= 90.0) / total * 100.0)
    cf = float(np.sum((scores >= 70.0) & (scores < 90.0)) / total * 100.0)
    lw = float(np.sum((scores >= 50.0) & (scores < 70.0)) / total * 100.0)
    vl = float(np.sum(scores < 50.0) / total * 100.0)

    # Domain-specific summaries
    dom_summary: Dict[str, Dict[str, float]] = {}
    for domain, group in df.groupby("domain"):
        grp_scores = group["plddt"].values
        dom_summary[str(domain)] = {
            "mean_plddt": round(float(np.mean(grp_scores)), 2),
            "median_plddt": round(float(np.median(grp_scores)), 2),
            "pct_disorder": round(float(np.sum(grp_scores < 50.0) / len(grp_scores) * 100.0), 1),
            "residue_count": len(grp_scores),
        }

    summary = ConfidenceSummary(
        total_residues=total,
        mean_plddt=round(float(np.mean(scores)), 2),
        median_plddt=round(float(np.median(scores)), 2),
        pct_very_high=round(vh, 1),
        pct_confident=round(cf, 1),
        pct_low=round(lw, 1),
        pct_very_low=round(vl, 1),
        domain_metrics=dom_summary,
    )

    logger.info(
        f"AlphaFold evaluation: {total} residues, Mean pLDDT: {summary.mean_plddt:.2f}, "
        f"Very High: {summary.pct_very_high}%, Disordered (<50): {summary.pct_very_low}%"
    )
    return df, summary
