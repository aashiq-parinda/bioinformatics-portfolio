import sys
from pathlib import Path

from typer.testing import CliRunner

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from src.cli import app

runner = CliRunner()


def test_cli_help() -> None:
    result = runner.invoke(app, ["--help"])
    assert result.exit_code == 0
    assert "Professional Protein Sequence & Homology Analysis" in result.stdout


def test_cli_analyze_command(tmp_path: Path) -> None:
    fasta_path = Path("projects/01-sequence-homology/examples/human_insulin.fasta")
    out_dir = tmp_path / "results"

    result = runner.invoke(app, ["analyze", str(fasta_path), "--output", str(out_dir)])
    assert result.exit_code == 0
    assert "Analyzing sequence from" in result.stdout
    assert "Physicochemical Properties" in result.stdout

    # Check generated files
    assert (out_dir / "sp_P01308_INS_HUMAN_statistics.json").exists() or any(
        out_dir.glob("*_statistics.json")
    )
    assert any(out_dir.glob("*.html"))
    assert (out_dir / "figures" / "amino_acid_composition.png").exists()
    assert (out_dir / "figures" / "hydropathy_profile.png").exists()
