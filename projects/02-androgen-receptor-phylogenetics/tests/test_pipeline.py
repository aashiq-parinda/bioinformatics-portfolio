"""End-to-end integration test for Project 02 pipeline."""

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
repo_root = project_root.parent.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from src.run_pipeline import run_pipeline


def test_pipeline_execution(tmp_path: Path) -> None:
    config_path = project_root / "configs" / "pipeline_config.json"
    result = run_pipeline(config_path)

    assert result["alignment"].exists()
    assert result["tree"].exists()
    assert result["html"].exists()
    assert result["json"].exists()
