"""Integration: Atlas API with real graph + embeddings on the compute server."""
from pathlib import Path
import pytest

from genome_atlas.api import Atlas


_HAS_ATLAS = Path.home().joinpath("pen-stack/data/graphs/atlas.gpickle").exists()


@pytest.mark.skipif(not _HAS_ATLAS, reason="Full atlas not available on this runner")
class TestAtlasFull:
    def test_query_system_returns_fields(self, demo_graph_path, demo_targets_path):
        atlas = Atlas(demo_graph_path, targets_path=demo_targets_path)
        result = atlas.query_system("SpCas9")
        assert "node_id" in result
        assert result.get("name") == "SpCas9"
        assert "mechanism_bucket" in result

    def test_query_system_all_14_systems(self, demo_graph_path, demo_targets_path):
        atlas = Atlas(demo_graph_path, targets_path=demo_targets_path)
        systems_df = atlas.systems()
        assert len(systems_df) == 14
        assert "name" in systems_df.columns

    def test_select_editor_top_k_shape(self, demo_graph_path, demo_targets_path):
        atlas = Atlas(demo_graph_path, targets_path=demo_targets_path)
        recs = atlas.select_editor("HEK293T", "deletion", 0, "AAV", True, top_k=5)
        assert 1 <= len(recs) <= 5
        for r in recs:
            assert 0.0 <= r.pen_score <= 1.0
            assert r.system != ""

    def test_similar_nodes_returns_dataframe(self, demo_graph_path,
                                             demo_embeddings_path, demo_targets_path):
        if not demo_embeddings_path.exists():
            pytest.skip("Embeddings not available")
        atlas = Atlas(demo_graph_path, demo_embeddings_path, demo_targets_path)
        sims = atlas.similar_nodes("System_1", node_type="System", top_k=3)
        assert len(sims) > 0
        assert "similarity" in sims.columns
        assert all(sims["similarity"] <= 1.0)
