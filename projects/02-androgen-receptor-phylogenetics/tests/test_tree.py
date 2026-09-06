"""Unit tests for phylogenetic tree reconstruction and Newick export."""

import sys
from pathlib import Path

project_root = Path(__file__).resolve().parent.parent
repo_root = project_root.parent.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(project_root) not in sys.path:
    sys.path.insert(0, str(project_root))

from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from src.tree import construct_phylogenetic_tree, export_newick_tree


def test_construct_phylogenetic_tree(tmp_path: Path) -> None:
    # Construct minimal test alignment of 4 sequences
    seqs = [
        SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWY"), id="Human"),
        SeqRecord(Seq("ACDEFGHIKLMNPQRSTVWY"), id="Chimp"),
        SeqRecord(Seq("ACDEFGHIKLMNPQRSSVWY"), id="Mouse"),
        SeqRecord(Seq("ACDEFGHIKAMNPQRSSVWY"), id="Zebrafish"),
    ]
    alignment = MultipleSeqAlignment(seqs)

    tree, dm = construct_phylogenetic_tree(alignment, model="identity", method="nj")
    assert tree is not None
    assert len(tree.get_terminals()) == 4

    # Test export
    out_nwk = tmp_path / "test.nwk"
    nwk_str = export_newick_tree(tree, out_nwk)
    assert out_nwk.exists()
    assert "Human" in nwk_str
    assert "Chimp" in nwk_str
