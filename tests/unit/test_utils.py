"""Tests for genome_atlas.utils.size - covers all branches including graph traversal."""
from __future__ import annotations

import networkx as nx
import pytest

from genome_atlas.api import Atlas
from genome_atlas.utils.size import system_total_size_aa


def _atlas_with_size(protein_aa: int) -> Atlas:
    """Build an Atlas with one System->Protein HAS_PROTEIN edge + length_map."""
    G = nx.MultiDiGraph()
    G.add_node("System_SpCas9", node_type="System", name="SpCas9")
    G.add_node("Protein_Q99ZW2", node_type="Protein", accession="Q99ZW2")
    G.add_edge("System_SpCas9", "Protein_Q99ZW2", edge_type="HAS_PROTEIN")

    atlas = Atlas.__new__(Atlas)
    atlas._G = G
    atlas._length_map = {"Protein_Q99ZW2": protein_aa}
    atlas._protein_map = {}
    atlas._system_by_name = {}
    atlas._domain_by_accession = {}
    atlas._protein_by_accession = {}
    atlas._rna_by_name = {}
    atlas._embeddings = None
    return atlas


class TestSystemTotalSizeAa:

    def test_returns_zero_when_no_graph(self):
        atlas = Atlas.__new__(Atlas)
        atlas._G = None
        atlas._length_map = {}
        result = system_total_size_aa(atlas, "System_SpCas9")
        assert result == 0

    def test_sums_protein_lengths(self):
        atlas = _atlas_with_size(1368)
        result = system_total_size_aa(atlas, "System_SpCas9")
        assert result == 1368

    def test_zero_for_unknown_protein(self):
        """Protein node exists in graph but not in length_map -> contributes 0."""
        atlas = _atlas_with_size(500)
        atlas._length_map = {}   # wipe map
        result = system_total_size_aa(atlas, "System_SpCas9")
        assert result == 0

    def test_sum_across_multiple_proteins(self):
        """Multi-protein system sums all subunit lengths."""
        G = nx.MultiDiGraph()
        G.add_node("System_CAST",    node_type="System", name="CAST")
        G.add_node("Protein_A",      node_type="Protein", accession="A")
        G.add_node("Protein_B",      node_type="Protein", accession="B")
        G.add_edge("System_CAST", "Protein_A", edge_type="HAS_PROTEIN")
        G.add_edge("System_CAST", "Protein_B", edge_type="HAS_PROTEIN")

        atlas = Atlas.__new__(Atlas)
        atlas._G = G
        atlas._length_map = {"Protein_A": 300, "Protein_B": 450}
        atlas._protein_map = {}
        atlas._system_by_name = {}
        atlas._domain_by_accession = {}
        atlas._protein_by_accession = {}
        atlas._rna_by_name = {}
        atlas._embeddings = None

        result = system_total_size_aa(atlas, "System_CAST")
        assert result == 750

    def test_ignores_non_has_protein_edges(self):
        """HAS_RNA and HAS_DOMAIN edges should not contribute to size."""
        G = nx.MultiDiGraph()
        G.add_node("System_SpCas9", node_type="System", name="SpCas9")
        G.add_node("RNA_sgRNA",     node_type="RNA",    name="sgRNA")
        G.add_edge("System_SpCas9", "RNA_sgRNA", edge_type="HAS_RNA")

        atlas = Atlas.__new__(Atlas)
        atlas._G = G
        atlas._length_map = {"RNA_sgRNA": 999}  # would corrupt size if included
        atlas._protein_map = {}
        atlas._system_by_name = {}
        atlas._domain_by_accession = {}
        atlas._protein_by_accession = {}
        atlas._rna_by_name = {}
        atlas._embeddings = None

        result = system_total_size_aa(atlas, "System_SpCas9")
        assert result == 0
