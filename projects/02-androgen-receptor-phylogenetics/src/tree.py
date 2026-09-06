"""Phylogenetic tree reconstruction using Distance Matrix and Neighbor-Joining."""

import io
from pathlib import Path
from typing import Tuple

from Bio import Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import (
    DistanceCalculator,
    DistanceTreeConstructor,
    _DistanceMatrix,
)

from shared.logging.logger import get_logger

logger = get_logger("projects.02.tree")


def construct_phylogenetic_tree(
    alignment: MultipleSeqAlignment,
    model: str = "blosum62",
    method: str = "nj",
) -> Tuple[Phylo.BaseTree.Tree, _DistanceMatrix]:
    """Construct a phylogenetic tree from an alignment using Neighbor-Joining or UPGMA.

    Args:
        alignment: Biopython MultipleSeqAlignment object.
        model: Protein distance model ('blosum62', 'identity').
        method: Tree construction method ('nj' for Neighbor Joining, 'upgma' for UPGMA).

    Returns:
        Tuple of (Constructed Tree object, DistanceMatrix).
    """
    calculator = DistanceCalculator(model)
    dm = calculator.get_distance(alignment)

    constructor = DistanceTreeConstructor(calculator, method)
    tree = constructor.nj(dm) if method == "nj" else constructor.upgma(dm)

    # Root at the outgroup (Zebrafish or longest branch) if present
    taxa_names = [record.id for record in alignment]
    outgroups = [
        name for name in taxa_names if "Danio" in name or "DANRE" in name or "rerio" in name
    ]
    if outgroups:
        try:
            tree.root_with_outgroup({"name": outgroups[0]})
            logger.info(f"Rooted phylogenetic tree with outgroup: {outgroups[0]}")
        except Exception as e:
            logger.warning(f"Could not root with outgroup: {e}")

    logger.info(f"Constructed {method.upper()} tree for {len(taxa_names)} taxa.")
    return tree, dm


def export_newick_tree(tree: Phylo.BaseTree.Tree, output_path: str | Path) -> str:
    """Export phylogenetic tree in standard Newick format."""
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:
        Phylo.write(tree, f, "newick")

    # Also return the Newick string
    buf = io.StringIO()
    Phylo.write(tree, buf, "newick")
    return buf.getvalue().strip()
