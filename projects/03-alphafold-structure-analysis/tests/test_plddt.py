"""Unit tests for AlphaFold pLDDT extraction and confidence metric calculations."""

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
repo_root = project_root.parent.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from src.alphafold_metrics import extract_alphafold_plddt


def test_extract_alphafold_plddt() -> None:
    af_pdb = project_root / "data" / "af_p10275_human_ar.pdb"
    domains = {
        "NTD": {"name": "NTD", "start": 1, "end": 555},
        "LBD": {"name": "LBD", "start": 676, "end": 919},
    }

    df, summary = extract_alphafold_plddt(af_pdb, domain_definitions=domains)

    assert len(df) in (919, 920)
    assert "plddt" in df.columns
    assert "confidence_tier" in df.columns
    assert summary.total_residues in (919, 920)

    # NTD should have significant disorder (<50 pLDDT)
    ntd_metrics = summary.domain_metrics.get("NTD", {})
    assert ntd_metrics.get("pct_disorder", 0) > 50.0

    # LBD should have high confidence (>70 mean pLDDT)
    lbd_metrics = summary.domain_metrics.get("LBD", {})
    assert lbd_metrics.get("mean_plddt", 0) > 75.0
