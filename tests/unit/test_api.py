"""Unit tests for Atlas API with mock graph."""
import pickle
from pathlib import Path

import networkx as nx
import pytest

from genome_atlas.api import Atlas


def _make_mock_graph():
    G = nx.MultiDiGraph()
    G.add_node("System_SpCas9", node_type="System", name="SpCas9")
    G.add_node("Protein_sp|P00000|TEST_HUMAN", node_type="Protein", name="TEST")
    G.add_node("Domain_PF00000", node_type="Domain", name="PF00000")
    G.add_edge("System_SpCas9", "Protein_sp|P00000|TEST_HUMAN", edge_type="HAS_PROTEIN")
    G.add_edge("Protein_sp|P00000|TEST_HUMAN", "Domain_PF00000", edge_type="HAS_DOMAIN")
    return G


class TestAtlasMock:
    def test_query_system(self, tmp_path):
        gpath = tmp_path / "test.gpickle"
        with gpath.open("wb") as f:
            pickle.dump(_make_mock_graph(), f)
        atlas = Atlas(gpath)
        result = atlas.query_system("SpCas9")
        assert result["node_id"] == "System_SpCas9"
        assert result["name"] == "SpCas9"

    def test_query_system_missing(self, tmp_path):
        gpath = tmp_path / "test.gpickle"
        with gpath.open("wb") as f:
            pickle.dump(_make_mock_graph(), f)
        atlas = Atlas(gpath)
        with pytest.raises(KeyError):
            atlas.query_system("Missing")

    def test_domains_of_protein(self, tmp_path):
        gpath = tmp_path / "test.gpickle"
        with gpath.open("wb") as f:
            pickle.dump(_make_mock_graph(), f)
        atlas = Atlas(gpath)
        domains = atlas.domains_of_protein("sp|P00000|TEST_HUMAN")
        assert len(domains) == 1
        assert domains.iloc[0]["node_id"] == "Domain_PF00000"

    def test_systems_list(self, tmp_path):
        gpath = tmp_path / "test.gpickle"
        with gpath.open("wb") as f:
            pickle.dump(_make_mock_graph(), f)
        atlas = Atlas(gpath)
        systems = atlas.systems()
        assert len(systems) == 1
