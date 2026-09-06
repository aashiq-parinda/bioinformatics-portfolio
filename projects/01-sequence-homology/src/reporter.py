"""Automated report generation and visualization assembler for sequence analysis."""

from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

from shared.io.tabular import save_dataframe, save_json
from shared.logging.logger import get_logger
from shared.reporting.html_generator import render_html_report
from shared.visualization.plots import (
    plot_amino_acid_composition,
    plot_blast_hits_distribution,
    plot_hydropathy_profile,
)

from .blast import BlastHit, hits_to_dataframe
from .stats import (
    SequenceStatistics,
    compute_hydropathy_profile,
)

logger = get_logger("projects.01.reporter")


def generate_analysis_report(
    stats: SequenceStatistics,
    sequence: str,
    output_dir: str | Path,
    hits: Optional[List[BlastHit]] = None,
) -> Dict[str, Path]:
    """Generate all analytical outputs: JSON summary, CSV hits, figures, and HTML report.

    Args:
        stats: Computed sequence statistics.
        sequence: Full raw sequence string.
        output_dir: Target output directory.
        hits: Optional list of filtered BLAST hits.

    Returns:
        Dictionary mapping artifact names to their generated file paths.
    """
    out_dir = Path(output_dir)
    figs_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    figs_dir.mkdir(parents=True, exist_ok=True)

    generated_files: Dict[str, Path] = {}

    from shared.utils.validators import sanitize_filename

    safe_id = sanitize_filename(stats.sequence_id)

    # 1. Save JSON statistics
    json_path = out_dir / f"{safe_id}_statistics.json"
    full_data = {
        "metadata": {
            "analysis_date": datetime.now().isoformat(),
            "pipeline": "bioseq-homology-v1.0",
        },
        "sequence_statistics": stats.to_dict(),
        "blast_hits_count": len(hits) if hits else 0,
    }
    save_json(full_data, json_path)
    generated_files["json"] = json_path

    # 2. Visualizations
    comp_fig = figs_dir / "amino_acid_composition.png"
    plot_amino_acid_composition(
        stats.amino_acid_composition_pct,
        comp_fig,
        title=f"Amino Acid Composition ({stats.sequence_id})",
    )
    generated_files["composition_plot"] = comp_fig

    hydropathy_scores = compute_hydropathy_profile(sequence, window_size=9)
    hydro_fig = figs_dir / "hydropathy_profile.png"
    plot_hydropathy_profile(
        hydropathy_scores,
        window_size=9,
        output_path=hydro_fig,
        title=f"Kyte-Doolittle Hydropathy Profile ({stats.sequence_id})",
    )
    generated_files["hydropathy_plot"] = hydro_fig

    # 3. BLAST Table & Distribution Figure
    hits_df = hits_to_dataframe(hits or [])
    csv_path = out_dir / f"{safe_id}_blast_hits.csv"
    save_dataframe(hits_df, csv_path, sep=",")
    generated_files["csv"] = csv_path

    blast_fig = figs_dir / "blast_hits_distribution.png"
    plot_blast_hits_distribution(
        hits_df,
        blast_fig,
        title=f"BLAST Homology Hit Distribution ({stats.sequence_id})",
    )
    generated_files["blast_plot"] = blast_fig

    # 4. Build HTML Sections and Stats Cards
    stat_cards = [
        {"label": "Sequence Length", "value": f"{stats.length:,} aa"},
        {"label": "Molecular Weight", "value": f"{stats.molecular_weight_da / 1000:.1f} kDa"},
        {"label": "Isoelectric Point (pI)", "value": f"{stats.isoelectric_point:.2f}"},
        {"label": "GRAVY Hydropathy", "value": f"{stats.gravy_hydropathy:.3f}"},
        {"label": "Aromaticity", "value": f"{stats.aromaticity:.3f}"},
        {
            "label": "Instability Index",
            "value": f"{stats.instability_index:.1f} ({'Stable' if stats.is_stable else 'Unstable'})",
        },
    ]

    # Format Top BLAST hits table for HTML report
    table_cols = ["Accession", "Organism", "Identity %", "Query Cov %", "E-value", "Bit Score"]
    table_rows = []
    if hits:
        for h in hits[:10]:
            table_rows.append(
                [
                    h.accession,
                    h.organism,
                    f"{h.identity_pct:.1f}%",
                    f"{h.query_coverage_pct:.1f}%",
                    f"{h.evalue:.2e}" if h.evalue < 0.001 else f"{h.evalue:.3f}",
                    f"{h.bit_score:.1f}",
                ]
            )

    sections = [
        {
            "title": "1. Physicochemical Properties & Hydropathy",
            "content": (
                f"Protein identifier <strong>{stats.sequence_id}</strong>: {stats.description}. "
                f"The sequence possesses a computed molecular weight of {stats.molecular_weight_da:.2f} Da and "
                f"a theoretical isoelectric point of {stats.isoelectric_point:.2f}. "
                f"The GRAVY index of {stats.gravy_hydropathy} characterizes the overall hydropathic nature."
            ),
            "figure_url": "figures/hydropathy_profile.png",
            "figure_caption": "Sliding window hydropathy profile (Kyte-Doolittle scale, window size = 9).",
        },
        {
            "title": "2. Amino Acid Residue Composition",
            "content": "Relative frequency distribution of standard amino acid residues.",
            "figure_url": "figures/amino_acid_composition.png",
            "figure_caption": "Per-residue percentage frequency distribution.",
        },
        {
            "title": "3. Homology & BLAST Analysis Summary",
            "content": f"A total of {len(hits) if hits else 0} homologous sequences met all identity and E-value significance cutoffs.",
            "figure_url": "figures/blast_hits_distribution.png" if hits else None,
            "figure_caption": "Identity % vs. -log10(E-value) distribution of top homologous sequences.",
            "table_columns": table_cols if hits else None,
            "table_data": table_rows if hits else None,
        },
    ]

    html_path = out_dir / f"{safe_id}_analysis_report.html"
    render_html_report(
        title=f"Sequence & Homology Report: {stats.sequence_id}",
        subtitle=stats.description or "Automated Bioinformatic Sequence Assessment",
        pipeline_name="bioseq analyze & blast",
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        stats=stat_cards,
        sections=sections,
        disclaimer="Computational sequence analysis is provided for scientific research and hypothesis generation only.",
        output_path=html_path,
    )
    generated_files["html"] = html_path

    logger.info(f"Report generated successfully at {html_path}")
    return generated_files
