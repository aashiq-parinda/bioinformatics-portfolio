"""Command Line Interface (CLI) for Automated Protein Sequence & Homology Analysis."""

import io
import sys
from pathlib import Path
from typing import Optional

import typer
from rich.console import Console
from rich.table import Table

# Allow running directly as a script
if __name__ == "__main__" and __package__ is None:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from src.blast import (
        filter_blast_hits,
        parse_blast_xml_stream,
        run_local_blast,
        run_remote_ncbi_blast,
    )
    from src.reporter import generate_analysis_report
    from src.stats import calculate_sequence_statistics
else:
    try:
        from .blast import (
            filter_blast_hits,
            parse_blast_xml_stream,
            run_local_blast,
            run_remote_ncbi_blast,
        )
        from .reporter import generate_analysis_report
        from .stats import calculate_sequence_statistics
    except ImportError:
        from src.blast import (
            filter_blast_hits,
            parse_blast_xml_stream,
            run_local_blast,
            run_remote_ncbi_blast,
        )
        from src.reporter import generate_analysis_report
        from src.stats import calculate_sequence_statistics
from shared.io.fasta import parse_fasta
from shared.logging.logger import get_logger

app = typer.Typer(
    name="bioseq",
    help="Professional Protein Sequence & Homology Analysis CLI Tool",
    add_completion=False,
)
console = Console()
logger = get_logger("bioseq.cli")


@app.command("analyze")
def analyze(
    fasta_path: Path = typer.Argument(..., help="Path to input protein FASTA file"),
    output: Path = typer.Option(
        Path("./results"), "--output", "-o", help="Output directory for results"
    ),
) -> None:
    """Perform comprehensive physicochemical sequence analysis on input FASTA."""
    console.print(f"[bold cyan]🧬 Analyzing sequence from:[/bold cyan] {fasta_path}")
    entries = parse_fasta(fasta_path)
    output.mkdir(parents=True, exist_ok=True)

    for entry in entries:
        console.print(
            f"\n[bold green]Processing record:[/bold green] {entry.id} ({entry.length} aa)"
        )
        stats = calculate_sequence_statistics(
            entry.sequence, seq_id=entry.id, description=entry.description
        )

        # Print Rich Summary Table to Console
        table = Table(title=f"Physicochemical Properties: {entry.id}")
        table.add_column("Property", style="cyan", justify="left")
        table.add_column("Value", style="magenta", justify="right")
        table.add_row("Sequence Length", f"{stats.length} aa")
        table.add_row("Molecular Weight", f"{stats.molecular_weight_da:,.2f} Da")
        table.add_row("Theoretical pI", f"{stats.isoelectric_point:.2f}")
        table.add_row("GRAVY Hydropathy", f"{stats.gravy_hydropathy:+.3f}")
        table.add_row("Aromaticity", f"{stats.aromaticity:.3f}")
        table.add_row(
            "Instability Index",
            f"{stats.instability_index:.2f} ({'Stable' if stats.is_stable else 'Unstable'})",
        )
        console.print(table)

        # Generate Full Reports
        files = generate_analysis_report(stats, entry.sequence, output_dir=output)
        console.print(f"[bold green]✓ Artifacts generated in:[/bold green] {output.resolve()}")
        console.print(f"  • HTML Report: [link={files['html']}]{files['html'].name}[/link]")
        console.print(f"  • JSON Statistics: {files['json'].name}")


@app.command("blast")
def blast(
    fasta_path: Path = typer.Argument(..., help="Path to input protein FASTA file"),
    database: str = typer.Option(
        "swissprot", "--database", "-db", help="BLAST database (e.g., swissprot, nr)"
    ),
    evalue: float = typer.Option(1e-5, "--evalue", "-e", help="Maximum E-value threshold"),
    min_identity: float = typer.Option(
        30.0, "--min-identity", "-id", help="Minimum percent identity (0-100)"
    ),
    min_coverage: float = typer.Option(
        40.0, "--min-coverage", "-cov", help="Minimum query coverage % (0-100)"
    ),
    max_hits: int = typer.Option(20, "--max-hits", "-n", help="Maximum hits to retain"),
    local_db: Optional[str] = typer.Option(
        None, "--local-db", help="Path to local BLAST database (if using local blastp)"
    ),
    output: Path = typer.Option(
        Path("./results"), "--output", "-o", help="Output directory for results"
    ),
) -> None:
    """Execute BLAST homology search (remote NCBI QBLAST or local BLAST+), filter hits, and generate reports."""
    console.print(f"[bold cyan]🔍 Executing BLAST for:[/bold cyan] {fasta_path}")
    entries = parse_fasta(fasta_path)
    output.mkdir(parents=True, exist_ok=True)

    for entry in entries:
        console.print(
            f"\n[bold yellow]Query Sequence:[/bold yellow] {entry.id} ({entry.length} aa)"
        )
        stats = calculate_sequence_statistics(
            entry.sequence, seq_id=entry.id, description=entry.description
        )

        if local_db:
            console.print(f"[blue]Running local blastp against {local_db}...[/blue]")
            xml_str = run_local_blast(
                fasta_path, local_db, evalue=evalue, max_target_seqs=max_hits * 2
            )
        else:
            console.print(
                f"[blue]Querying NCBI QBLAST (database={database}, evalue={evalue})...[/blue]"
            )
            xml_str = run_remote_ncbi_blast(
                entry.sequence,
                database=database,
                evalue_threshold=evalue,
                hitlist_size=max_hits * 2,
            )

        raw_hits = parse_blast_xml_stream(io.StringIO(xml_str), query_length=entry.length)
        filtered_hits = filter_blast_hits(
            raw_hits,
            min_identity=min_identity,
            max_evalue=evalue,
            min_query_coverage=min_coverage,
            max_hits=max_hits,
        )

        console.print(
            f"[green]Identified {len(raw_hits)} raw hits; {len(filtered_hits)} passed filtering criteria.[/green]"
        )

        # Print Top Hits Table
        hit_table = Table(title=f"Top Homologs: {entry.id}")
        hit_table.add_column("Accession", style="cyan")
        hit_table.add_column("Organism", style="white")
        hit_table.add_column("Identity %", style="green", justify="right")
        hit_table.add_column("Query Cov %", style="yellow", justify="right")
        hit_table.add_column("E-value", style="magenta", justify="right")
        hit_table.add_column("Bit Score", style="blue", justify="right")

        for h in filtered_hits[:10]:
            hit_table.add_row(
                h.accession,
                h.organism[:25],
                f"{h.identity_pct:.1f}%",
                f"{h.query_coverage_pct:.1f}%",
                f"{h.evalue:.2e}" if h.evalue < 0.001 else f"{h.evalue:.3f}",
                f"{h.bit_score:.1f}",
            )
        console.print(hit_table)

        files = generate_analysis_report(
            stats, entry.sequence, output_dir=output, hits=filtered_hits
        )
        console.print(
            f"[bold green]✓ Complete analysis generated in:[/bold green] {output.resolve()}"
        )
        console.print(f"  • HTML Report: {files['html'].name}")


def main() -> None:
    app()


if __name__ == "__main__":
    main()
