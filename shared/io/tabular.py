"""Tabular data input/output utilities with Pandas, JSON, and TSV/CSV helpers."""

import json
from pathlib import Path
from typing import Any, Dict, List

import pandas as pd

from shared.logging.logger import get_logger

logger = get_logger("shared.io.tabular")


def load_dataframe(file_path: str | Path, sep: str = "\t") -> pd.DataFrame:
    """Load a TSV/CSV tabular file into a pandas DataFrame."""
    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"Tabular file not found: {path}")

    # Detect delimiter automatically if not specified
    if path.suffix == ".csv":
        sep = ","
    df = pd.read_csv(path, sep=sep)
    logger.info(f"Loaded DataFrame with shape {df.shape} from {path.name}")
    return df


def save_dataframe(
    df: pd.DataFrame, output_path: str | Path, sep: str = "\t", index: bool = False
) -> None:
    """Save a DataFrame to TSV or CSV format."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".csv":
        sep = ","
    df.to_csv(path, sep=sep, index=index)
    logger.info(f"Saved DataFrame to {path}")


def save_json(data: Dict[str, Any] | List[Any], output_path: str | Path, indent: int = 2) -> None:
    """Save dictionary or list data to a formatted JSON file."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=indent, default=str)
    logger.info(f"Saved JSON data to {path}")


def load_json(file_path: str | Path) -> Any:
    """Load data from a JSON file."""
    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"JSON file not found: {path}")
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)
