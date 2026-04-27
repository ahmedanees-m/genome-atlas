"""Integration: Atlas API with graph + embeddings."""
from pathlib import Path
import pytest

from genome_atlas.api import Atlas


@pytest.mark.skipif(
    not Path.home().joinpath("pen-stack/data/graphs/atlas.gpickle").exists(),
    reason="Full atlas not available on this runner",
)
class TestAtlasFull:
    def test_query_system_returns_fields(self, demo_graph_path, demo_targets_path):
        atlas = Atlas(str(demo_graph_path), targets_path=str(demo_targets_path))
        result = atlas.query_system("System_SpCas9")
        assert result["node_id"] == "System_SpCas9"
        assert "name" in result
        assert "domains" in result or "mechanisms" in result

    def test_select_editor_top_k_shape(self, demo_graph_path, demo_targets_path):
        atlas = Atlas(str(demo_graph_path), targets_path=str(demo_targets_path))
        recs = atlas.select_editor("HEK293T", "deletion", 0, "AAV", True, top_k=5)
        assert 1 <= len(recs) <= 5
        for r in recs:
            assert 0.0 <= r.pen_score <= 1.0

    def test_similar_nodes_nonempty(self, demo_graph_path, demo_embeddings_path, demo_targets_path):
        if not demo_embeddings_path.exists():
            pytest.skip("Embeddings not available")
        atlas = Atlas(str(demo_graph_path), str(demo_embeddings_path), str(demo_targets_path))
        sims = atlas.similar_nodes("Protein_sp|P10178|NKX31_HUMAN", top_k=5)
        assert len(sims) > 0
        assert all("pen_score" in s for s in sims)
