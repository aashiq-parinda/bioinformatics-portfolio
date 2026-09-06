"""Over-Representation Analysis (ORA) for Gene Ontology and pathway gene sets."""

from typing import Any, Dict, List, Set

import numpy as np
import pandas as pd
from scipy import stats

from shared.logging.logger import get_logger

logger = get_logger("projects.05.pathways")

# Curated Biological Process and Pathway Gene Sets (MSigDB Hallmark / GO BP / KEGG)
REFERENCE_GENE_SETS: Dict[str, Dict[str, Any]] = {
    "GO:0048545": {
        "term_id": "GO:0048545",
        "name": "Response to Steroid Hormone",
        "category": "Gene Ontology BP",
        "genes": {
            "KLK3",
            "FKBP5",
            "TMPRSS2",
            "NKX3-1",
            "SLC45A3",
            "PMEPA1",
            "IGF1R",
            "MYC",
            "BCL2",
            "FOXA1",
        },
    },
    "KEGG:05215": {
        "term_id": "KEGG:05215",
        "name": "Androgen Receptor Signaling Pathway",
        "category": "KEGG Pathway",
        "genes": {
            "KLK3",
            "FKBP5",
            "TMPRSS2",
            "NKX3-1",
            "SLC45A3",
            "ELL2",
            "ABCC4",
            "STEAP1",
            "STEAP2",
            "NDRG1",
        },
    },
    "GO:0006629": {
        "term_id": "GO:0006629",
        "name": "Lipid Metabolic Process",
        "category": "Gene Ontology BP",
        "genes": {"ACSL3", "PPAP2A", "SEC14L2", "GUCY1A3", "MAOA", "STEAP1"},
    },
    "GO:0008283": {
        "term_id": "GO:0008283",
        "name": "Cell Proliferation & Growth Control",
        "category": "Gene Ontology BP",
        "genes": {"MYC", "MKI67", "CDKN1A", "CCND1", "AURKA", "TOP2A", "IGF1R", "TGFB2"},
    },
    "HALLMARK:01": {
        "term_id": "HALLMARK:01",
        "name": "Androgen Response Hallmark",
        "category": "MSigDB Hallmark",
        "genes": {
            "KLK3",
            "FKBP5",
            "TMPRSS2",
            "NKX3-1",
            "SLC45A3",
            "PMEPA1",
            "STEAP1",
            "NDRG1",
            "ACSL3",
        },
    },
    "GO:0006915": {
        "term_id": "GO:0006915",
        "name": "Regulation of Apoptosis",
        "category": "Gene Ontology BP",
        "genes": {"BCL2", "CDKN1A", "MYC", "IGFBP3", "IL6"},
    },
}


def run_pathway_enrichment(
    significant_genes: Set[str],
    background_genes: Set[str],
    p_val_cutoff: float = 0.05,
) -> pd.DataFrame:
    """Perform Over-Representation Analysis using Fisher's exact hypergeometric test.

    Args:
        significant_genes: Set of gene IDs passing significance thresholds.
        background_genes: Total background universe of tested genes.
        p_val_cutoff: Raw p-value significance threshold.

    Returns:
        DataFrame of enriched pathways ranked by significance.
    """
    total_bg = len(background_genes)  # Total background population
    n_sig = len(significant_genes & background_genes)  # Number of significant genes in universe

    if n_sig == 0 or total_bg == 0:
        return pd.DataFrame()

    results: List[Dict[str, Any]] = []

    for term_id, term_info in REFERENCE_GENE_SETS.items():
        set_genes = term_info["genes"] & background_genes
        pathway_pop = len(set_genes)  # Total genes in the pathway

        if pathway_pop < 2:
            continue

        # Overlap between significant genes and pathway
        overlap_genes = significant_genes & set_genes
        k = len(overlap_genes)

        if k == 0:
            continue

        # Hypergeometric test: probability of drawing >= k successes
        p_value = float(stats.hypergeom.sf(k - 1, total_bg, pathway_pop, n_sig))

        # Fold enrichment = (k / n_sig) / (pathway_pop / total_bg)
        expected_overlap = (n_sig * pathway_pop) / total_bg
        fold_enrichment = (k / expected_overlap) if expected_overlap > 0 else 0.0

        if p_value <= p_val_cutoff:
            results.append(
                {
                    "term_id": term_id,
                    "pathway_name": term_info["name"],
                    "category": term_info["category"],
                    "overlap_count": k,
                    "pathway_size": pathway_pop,
                    "fold_enrichment": round(fold_enrichment, 2),
                    "p_value": p_value,
                    "neg_log10_p": round(-np.log10(p_value + 1e-300), 2),
                    "overlapping_genes": ", ".join(sorted(overlap_genes)),
                }
            )

    df = pd.DataFrame(results)
    if not df.empty:
        df = df.sort_values("p_value").reset_index(drop=True)

    logger.info(f"Pathway enrichment identified {len(df)} significant terms at p < {p_val_cutoff}.")
    return df
