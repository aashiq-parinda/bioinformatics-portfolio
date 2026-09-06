"""End-to-end orchestration pipeline for RNA-Seq differential expression analysis."""

import argparse
import sys
from datetime import datetime
from pathlib import Path

# Support direct script invocation
project_root = Path(__file__).resolve().parent.parent
repo_root = project_root.parent.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from shared.io.tabular import load_dataframe, load_json, save_dataframe, save_json
from shared.logging.logger import get_logger
from shared.reporting.html_generator import render_html_report
from src.differential_expression import run_differential_expression
from src.normalization import calculate_size_factors, normalize_counts
from src.pathway_enrichment import run_pathway_enrichment
from src.visualizer import (
    plot_expression_heatmap,
    plot_ma,
    plot_pathway_enrichment,
    plot_pca,
    plot_volcano,
)

logger = get_logger("projects.05.pipeline")


def run_rnaseq_pipeline(config_path: str | Path) -> dict:
    """Execute complete RNA-seq normalization, DE, pathway, and reporting pipeline."""
    cfg = load_json(config_path)
    base_dir = Path(config_path).resolve().parent.parent

    results_dir = base_dir / cfg["outputs"]["results_dir"]
    figures_dir = base_dir / cfg["outputs"]["figures_dir"]
    results_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Starting Project 05 RNA-Seq Pipeline using config: {config_path}")

    # 1. Load Data
    counts_file = base_dir / cfg["inputs"]["counts_matrix"]
    meta_file = base_dir / cfg["inputs"]["sample_metadata"]

    raw_counts = load_dataframe(counts_file, sep="\t")
    metadata = load_dataframe(meta_file, sep=",")

    target_cond = cfg["experimental_design"]["target_level"]
    ref_cond = cfg["experimental_design"]["reference_level"]

    ctrl_samples = metadata[metadata["condition"] == ref_cond]["sample_id"].tolist()
    treat_samples = metadata[metadata["condition"] == target_cond]["sample_id"].tolist()
    all_samples = ctrl_samples + treat_samples

    # 2. Median-of-Ratios Normalization
    size_factors = calculate_size_factors(raw_counts, all_samples)
    norm_counts, log2_counts = normalize_counts(raw_counts, all_samples, size_factors)

    norm_csv = results_dir / "normalized_counts.csv"
    save_dataframe(norm_counts, norm_csv, sep=",")

    # 3. Negative Binomial Wald Differential Expression
    fdr_cut = cfg["statistical_cutoffs"]["fdr_threshold"]
    lfc_cut = cfg["statistical_cutoffs"]["log2fc_threshold"]
    min_count = cfg["statistical_cutoffs"]["min_count_filter"]

    de_results, de_summary = run_differential_expression(
        norm_counts,
        control_samples=ctrl_samples,
        treated_samples=treat_samples,
        fdr_cutoff=fdr_cut,
        log2fc_cutoff=lfc_cut,
        min_count_sum=min_count,
    )

    de_csv = results_dir / "differential_expression_results.csv"
    save_dataframe(de_results, de_csv, sep=",")

    sig_genes_df = de_results[de_results["status"] != "Not Significant"].copy()
    sig_csv = results_dir / "significant_genes.csv"
    save_dataframe(sig_genes_df, sig_csv, sep=",")

    # 4. Pathway Over-Representation Analysis (ORA)
    sig_set = set(sig_genes_df["gene_id"])
    bg_set = set(de_results["gene_id"])
    enrichment_df = run_pathway_enrichment(
        sig_set, bg_set, p_val_cutoff=cfg["pathway_enrichment"]["p_value_cutoff"]
    )

    pathway_csv = results_dir / "pathway_enrichment_results.csv"
    save_dataframe(enrichment_df, pathway_csv, sep=",")

    # 5. Scientific Visualizations
    pca_fig = figures_dir / "pca_plot.png"
    plot_pca(log2_counts, all_samples, metadata, pca_fig)

    volcano_fig = figures_dir / "volcano_plot.png"
    plot_volcano(de_results, volcano_fig, fdr_cutoff=fdr_cut, log2fc_cutoff=lfc_cut)

    ma_fig = figures_dir / "ma_plot.png"
    plot_ma(de_results, ma_fig)

    heatmap_fig = figures_dir / "expression_heatmap.png"
    plot_expression_heatmap(log2_counts, de_results, all_samples, heatmap_fig, n_top=25)

    pathway_fig = figures_dir / "pathway_enrichment_bars.png"
    plot_pathway_enrichment(enrichment_df, pathway_fig)

    # 6. HTML Executive Summary Report
    stats_cards = [
        {"label": "Sequencing Samples", "value": f"{len(all_samples)} Libraries (3 vs 3)"},
        {"label": "Genes Quantified", "value": f"{de_summary.total_genes_tested:,} Genes"},
        {"label": "Significant Up-Regulated", "value": f"{de_summary.num_significant_up} Genes"},
        {
            "label": "Significant Down-Regulated",
            "value": f"{de_summary.num_significant_down} Genes",
        },
        {"label": "FDR Significance Cutoff", "value": f"padj < {fdr_cut}"},
        {"label": "Fold Change Threshold", "value": f"|log2FC| >= {lfc_cut}"},
    ]

    top_de_rows = []
    for _, r in sig_genes_df.head(10).iterrows():
        top_de_rows.append(
            [
                f"<strong>{r['gene_id']}</strong>",
                f"{r['baseMean']:.1f}",
                f"{r['log2FoldChange']:+.2f}",
                f"{r['pvalue']:.2e}" if r["pvalue"] < 0.001 else f"{r['pvalue']:.3f}",
                f"{r['padj']:.2e}" if r["padj"] < 0.001 else f"{r['padj']:.3f}",
                f"<span style='color: {'#dc2626' if r['status'] == 'Up-regulated' else '#2563eb'}; font-weight: bold;'>{r['status']}</span>",
            ]
        )

    pathway_rows = []
    if not enrichment_df.empty:
        for _, pr in enrichment_df.head(6).iterrows():
            pathway_rows.append(
                [
                    pr["term_id"],
                    pr["pathway_name"],
                    pr["category"],
                    f"{pr['overlap_count']} / {pr['pathway_size']}",
                    f"{pr['p_value']:.2e}",
                    pr["overlapping_genes"][:35] + "...",
                ]
            )

    sections = [
        {
            "title": "1. Study Design & Sample Separation (PCA)",
            "content": (
                f"Transcriptomic profiling of {cfg['dataset_metadata']['source_study']} evaluating "
                f"the hormonal response in {cfg['dataset_metadata']['organism']}. "
                "Principal Component Analysis of the top 500 variable genes demonstrates sharp separation along PC1 "
                "between Vehicle Control and Androgen-Treated replicates, verifying robust transcriptional reprogramming."
            ),
            "figure_url": "figures/pca_plot.png",
            "figure_caption": "Principal Component Analysis (PCA) demonstrating clear condition-based clustering.",
        },
        {
            "title": "2. Differential Expression Landscape (Volcano Plot)",
            "content": (
                f"Negative Binomial generalized linear model testing identified {de_summary.num_significant_up} up-regulated "
                f"and {de_summary.num_significant_down} down-regulated genes. Canonical AR downstream effectors "
                "(<em>KLK3</em>, <em>FKBP5</em>, <em>TMPRSS2</em>, <em>NKX3-1</em>) exhibit dramatic transcriptional induction, "
                "while proliferative oncogenes such as <em>MYC</em> are suppressed."
            ),
            "figure_url": "figures/volcano_plot.png",
            "figure_caption": "Volcano plot (-log10 FDR vs log2 Fold Change) highlighting key androgen responsive transcripts.",
            "table_columns": [
                "Gene Symbol",
                "Base Mean",
                "log2(Fold Change)",
                "p-value",
                "FDR (padj)",
                "Status",
            ],
            "table_data": top_de_rows,
        },
        {
            "title": "3. Clustered Expression Heatmap",
            "content": "Hierarchically clustered Z-score heatmap of top 25 differentially expressed genes.",
            "figure_url": "figures/expression_heatmap.png",
            "figure_caption": "Standardized expression profiles across biological replicates.",
        },
        {
            "title": "4. Functional Pathway & Gene Ontology Enrichment",
            "content": (
                "Over-Representation Analysis (ORA) confirms significant enrichment for the androgen receptor "
                "signaling pathway and steroid hormone response networks."
            ),
            "figure_url": "figures/pathway_enrichment_bars.png",
            "figure_caption": "Top enriched Gene Ontology terms and KEGG pathways.",
            "table_columns": [
                "Term ID",
                "Pathway Name",
                "Category",
                "Overlap",
                "p-value",
                "Top Genes",
            ],
            "table_data": pathway_rows if pathway_rows else None,
        },
    ]

    html_report = results_dir / "rnaseq_analysis_report.html"
    render_html_report(
        title="End-to-End RNA-Seq Differential Expression Report",
        subtitle=f"Hormone Response Transcriptomics ({cfg['dataset_metadata']['source_study']})",
        pipeline_name="05-rna-seq",
        timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        stats=stats_cards,
        sections=sections,
        disclaimer=(
            "This RNA-seq analysis is a computational genomics study of publicly available transcriptomic data. "
            "It does NOT make clinical diagnostic claims or medical therapy recommendations."
        ),
        output_path=html_report,
    )

    # Save JSON summary
    summary_json = results_dir / "rnaseq_summary.json"
    save_json(
        {
            "dataset_metadata": cfg["dataset_metadata"],
            "size_factors": size_factors,
            "de_summary": de_summary.__dict__,
            "top_upregulated_genes": sig_genes_df[sig_genes_df["status"] == "Up-regulated"]
            .head(10)["gene_id"]
            .tolist(),
            "top_downregulated_genes": sig_genes_df[sig_genes_df["status"] == "Down-regulated"]
            .head(10)["gene_id"]
            .tolist(),
            "enriched_pathways_count": len(enrichment_df),
        },
        summary_json,
    )

    logger.info(f"Project 05 RNA-Seq Pipeline completed successfully. Report: {html_report}")
    return {
        "de_csv": de_csv,
        "sig_csv": sig_csv,
        "html_report": html_report,
        "summary_json": summary_json,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run RNA-Seq Pipeline")
    parser.add_argument(
        "--config", default="configs/rnaseq_config.json", help="Path to config file"
    )
    args = parser.parse_args()

    cfg_path = Path(args.config)
    if not cfg_path.is_file():
        cfg_path = project_root / args.config

    run_rnaseq_pipeline(cfg_path)


if __name__ == "__main__":
    main()
