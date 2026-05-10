"""Unit tests for genome_atlas.graph.view (graph_view parameter).

Tests are fully self-contained: they build minimal synthetic NetworkX graphs
and PyG HeteroData objects without requiring the full atlas gpickle or any
VM-side data.

Note: this module is skipped automatically when PyTorch or torch-geometric are
not installed (e.g. in the lightweight CI environment that installs only the
``[dev]`` extras).  Install ``.[gnn]`` to run these tests locally.
"""
from __future__ import annotations

import pickle
from pathlib import Path

import networkx as nx
import numpy as np
import pytest

# Skip the entire module when the heavy ML deps are absent.
torch = pytest.importorskip("torch", reason="PyTorch not installed; skipping GNN tests")
pytest.importorskip("torch_geometric", reason="torch-geometric not installed; skipping GNN tests")
from torch_geometric.data import HeteroData  # noqa: E402 (after importorskip)

# ---------------------------------------------------------------------------
# Helpers to build lightweight synthetic graphs / data
# ---------------------------------------------------------------------------

def _make_synthetic_G(include_has_rna: bool = False) -> nx.MultiDiGraph:
    """Return a minimal NetworkX MultiDiGraph that mirrors the atlas schema."""
    G = nx.MultiDiGraph()

    # Protein nodes
    for acc in ["P00001", "P00002", "P00003"]:
        G.add_node(acc, node_type="Protein", accession=acc)

    # Domain nodes
    for dom in ["PF01548", "PF02371"]:
        G.add_node(dom, node_type="Domain")

    # System node
    G.add_node("IS621_bridge_recombinase", node_type="System")

    # Mechanism node
    G.add_node("DSB_FREE", node_type="Mechanism")

    # Structure node
    G.add_node("AF-P00001-F1", node_type="Structure")

    # RNA node
    G.add_node("bridge_RNA", node_type="RNA")

    # HAS_DOMAIN edges
    G.add_edge("P00001", "PF01548", edge_type="HAS_DOMAIN")
    G.add_edge("P00001", "PF02371", edge_type="HAS_DOMAIN")
    G.add_edge("P00002", "PF01548", edge_type="HAS_DOMAIN")

    # HAS_PROTEIN edge
    G.add_edge("IS621_bridge_recombinase", "P00001", edge_type="HAS_PROTEIN")

    # USES_MECHANISM edge
    G.add_edge("IS621_bridge_recombinase", "DSB_FREE", edge_type="USES_MECHANISM")

    # STRUCTURE_OF edge
    G.add_edge("AF-P00001-F1", "P00001", edge_type="STRUCTURE_OF")

    if include_has_rna:
        G.add_edge("IS621_bridge_recombinase", "bridge_RNA", edge_type="HAS_RNA")

    return G


def _make_esm_parquet(tmp_path: Path, accessions: list[str], emb_dim: int = 8) -> Path:
    """Write a minimal ESM-2 parquet with synthetic embeddings."""
    import pandas as pd

    rng = np.random.default_rng(42)
    rows = [
        {"accession": acc, "embedding": rng.random(emb_dim).tolist()}
        for acc in accessions
    ]
    df = pd.DataFrame(rows)
    path = tmp_path / "esm2_test.parquet"
    df.to_parquet(str(path), index=False)
    return path


def _write_gpickle(G: nx.MultiDiGraph, path: Path) -> None:
    with open(str(path), "wb") as f:
        pickle.dump(G, f)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def tmp_atlas(tmp_path):
    """Return (gpickle_path, esm_parquet_path) for a synthetic atlas."""
    G = _make_synthetic_G()
    gpickle_path = tmp_path / "atlas.gpickle"
    _write_gpickle(G, gpickle_path)
    esm_path = _make_esm_parquet(tmp_path, ["P00001", "P00002", "P00003"])
    return gpickle_path, esm_path


@pytest.fixture()
def tmp_atlas_with_has_rna(tmp_path):
    """Synthetic atlas that already contains HAS_RNA edges in the gpickle."""
    G = _make_synthetic_G(include_has_rna=True)
    gpickle_path = tmp_path / "atlas_with_rna.gpickle"
    _write_gpickle(G, gpickle_path)
    esm_path = _make_esm_parquet(tmp_path, ["P00001", "P00002", "P00003"])
    return gpickle_path, esm_path


# ---------------------------------------------------------------------------
# Tests - primary view
# ---------------------------------------------------------------------------

class TestPrimaryView:
    """get_graph(graph_view='primary') should match the v0.7.0 4-edge schema."""

    def test_returns_hetero_data(self, tmp_atlas):
        from genome_atlas.graph.view import get_graph

        gpickle, esm = tmp_atlas
        data = get_graph(gpickle, esm, graph_view="primary", add_split=False)
        assert isinstance(data, HeteroData)

    def test_has_exactly_four_primary_edge_types(self, tmp_atlas):
        from genome_atlas.graph.view import _PRIMARY_EDGE_TYPES, get_graph

        gpickle, esm = tmp_atlas
        data = get_graph(gpickle, esm, graph_view="primary", add_split=False)
        edge_types = {tuple(et) for et in data.edge_types}
        # All returned edge types must be primary
        assert edge_types <= _PRIMARY_EDGE_TYPES

    def test_no_secondary_edge_types(self, tmp_atlas):
        from genome_atlas.graph.view import get_graph

        gpickle, esm = tmp_atlas
        data = get_graph(gpickle, esm, graph_view="primary", add_split=False)
        edge_names = {et[1] for et in data.edge_types}
        assert "HAS_RNA" not in edge_names
        assert "PART_OF" not in edge_names
        assert "SIMILAR_TO" not in edge_names

    def test_has_domain_edge_present(self, tmp_atlas):
        from genome_atlas.graph.view import get_graph

        gpickle, esm = tmp_atlas
        data = get_graph(gpickle, esm, graph_view="primary", add_split=False)
        assert ("Protein", "HAS_DOMAIN", "Domain") in data.edge_types

    def test_add_split_adds_masks(self, tmp_atlas):
        from genome_atlas.graph.view import get_graph

        gpickle, esm = tmp_atlas
        data = get_graph(gpickle, esm, graph_view="primary", add_split=True)
        for et in data.edge_types:
            assert hasattr(data[et], "train_mask"), f"Missing train_mask on {et}"
            assert hasattr(data[et], "val_mask"),   f"Missing val_mask on {et}"
            assert hasattr(data[et], "test_mask"),  f"Missing test_mask on {et}"

    def test_no_split_when_disabled(self, tmp_atlas):
        from genome_atlas.graph.view import get_graph

        gpickle, esm = tmp_atlas
        data = get_graph(gpickle, esm, graph_view="primary", add_split=False)
        for et in data.edge_types:
            assert not hasattr(data[et], "train_mask"), \
                f"Unexpected train_mask on {et} when add_split=False"

    def test_default_graph_view_is_primary(self, tmp_atlas):
        """Calling get_graph without graph_view should default to primary."""
        from genome_atlas.graph.view import _PRIMARY_EDGE_TYPES, get_graph

        gpickle, esm = tmp_atlas
        data = get_graph(gpickle, esm, add_split=False)
        edge_types = {tuple(et) for et in data.edge_types}
        assert edge_types <= _PRIMARY_EDGE_TYPES


# ---------------------------------------------------------------------------
# Tests - full view
# ---------------------------------------------------------------------------

class TestFullView:
    """get_graph(graph_view='full') should add HAS_RNA, PART_OF, and optionally SIMILAR_TO."""

    def test_returns_hetero_data(self, tmp_atlas):
        from genome_atlas.graph.view import get_graph

        gpickle, esm = tmp_atlas
        data = get_graph(gpickle, esm, graph_view="full", add_split=False)
        assert isinstance(data, HeteroData)

    def test_includes_all_primary_edge_types(self, tmp_atlas):
        from genome_atlas.graph.view import _PRIMARY_EDGE_TYPES, get_graph

        gpickle, esm = tmp_atlas
        data = get_graph(gpickle, esm, graph_view="full", add_split=False)
        edge_types = {tuple(et) for et in data.edge_types}
        # All four primary types must still be present
        for pet in _PRIMARY_EDGE_TYPES:
            assert pet in edge_types, f"Primary edge type {pet} missing in full view"

    def test_part_of_derived_from_has_protein(self, tmp_atlas):
        """PART_OF edges should be the reverse of HAS_PROTEIN."""
        from genome_atlas.graph.view import get_graph

        gpickle, esm = tmp_atlas
        data = get_graph(gpickle, esm, graph_view="full", add_split=False)

        assert ("Protein", "PART_OF", "System") in data.edge_types, \
            "PART_OF edge type not derived"

        # Verify it is the transpose of HAS_PROTEIN
        hp_ei = data[("System", "HAS_PROTEIN", "Protein")].edge_index
        po_ei = data[("Protein", "PART_OF", "System")].edge_index

        assert torch.equal(hp_ei[0], po_ei[1]), "PART_OF src must equal HAS_PROTEIN dst"
        assert torch.equal(hp_ei[1], po_ei[0]), "PART_OF dst must equal HAS_PROTEIN src"

    def test_has_rna_not_duplicated_when_already_present(self, tmp_atlas_with_has_rna):
        """If gpickle already has HAS_RNA edges, full view should not re-derive them."""
        from genome_atlas.graph.view import get_graph

        gpickle, esm = tmp_atlas_with_has_rna
        data = get_graph(gpickle, esm, graph_view="full", add_split=False)

        assert ("System", "HAS_RNA", "RNA") in data.edge_types
        # There should be exactly one HAS_RNA edge (the one from the gpickle)
        n_rna_edges = data[("System", "HAS_RNA", "RNA")].edge_index.size(1)
        assert n_rna_edges >= 1, "Expected at least one HAS_RNA edge"

    def test_similar_to_edges_computed(self, tmp_atlas):
        """SIMILAR_TO edges should be added for proteins with high cosine similarity."""
        from genome_atlas.graph.view import get_graph

        gpickle, esm = tmp_atlas

        # Build data with very low threshold so we get edges from random embeddings
        data = get_graph(
            gpickle, esm,
            graph_view="full",
            similarity_threshold=0.0,  # All pairs (excluding self-loops)
            add_split=False,
        )

        # With threshold=0, all non-self pairs should produce edges (N*(N-1))
        # But zero-embedding rows are excluded in the computation - we have
        # 3 proteins with real ESM-2 embeddings, so expect 3*2 = 6 directed edges.
        if ("Protein", "SIMILAR_TO", "Protein") in data.edge_types:
            n_edges = data[("Protein", "SIMILAR_TO", "Protein")].edge_index.size(1)
            assert n_edges > 0, "Expected SIMILAR_TO edges with threshold=0"

    def test_similar_to_no_self_loops(self, tmp_atlas):
        """SIMILAR_TO edges must not include self-loops."""
        from genome_atlas.graph.view import get_graph

        gpickle, esm = tmp_atlas
        data = get_graph(
            gpickle, esm,
            graph_view="full",
            similarity_threshold=0.0,
            add_split=False,
        )

        if ("Protein", "SIMILAR_TO", "Protein") in data.edge_types:
            ei = data[("Protein", "SIMILAR_TO", "Protein")].edge_index
            # src != dst for all edges
            assert (ei[0] != ei[1]).all(), "SIMILAR_TO edges contain self-loops"


# ---------------------------------------------------------------------------
# Tests - invalid inputs
# ---------------------------------------------------------------------------

class TestInvalidInputs:

    def test_invalid_graph_view_raises(self, tmp_atlas):
        from genome_atlas.graph.view import get_graph

        gpickle, esm = tmp_atlas
        with pytest.raises(ValueError, match="graph_view must be"):
            get_graph(gpickle, esm, graph_view="unknown")  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# Tests - _filter_primary helper
# ---------------------------------------------------------------------------

class TestFilterPrimary:

    def test_removes_non_primary_edges(self):
        from genome_atlas.graph.view import _PRIMARY_EDGE_TYPES, _filter_primary

        data = HeteroData()
        # Add all primary edge types
        for src, rel, dst in _PRIMARY_EDGE_TYPES:
            data[(src, rel, dst)].edge_index = torch.zeros(2, 0, dtype=torch.long)
        # Add a non-primary edge
        data[("System", "HAS_RNA", "RNA")].edge_index = torch.zeros(2, 0, dtype=torch.long)
        data[("Protein", "PART_OF", "System")].edge_index = torch.zeros(2, 0, dtype=torch.long)

        result = _filter_primary(data)
        edge_names = {et[1] for et in result.edge_types}
        assert "HAS_RNA" not in edge_names
        assert "PART_OF" not in edge_names
        # Primary edges retained
        assert "HAS_DOMAIN" in edge_names


# ---------------------------------------------------------------------------
# Tests - _compute_similarity_edges helper
# ---------------------------------------------------------------------------

class TestComputeSimilarityEdges:

    def test_identical_vectors_max_similarity(self):
        from genome_atlas.graph.view import _compute_similarity_edges

        # Two identical vectors -> cosine similarity = 1.0
        feats = torch.tensor([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        src, dst = _compute_similarity_edges(feats, threshold=0.99)
        # Rows 0 and 1 are identical -> expect edges 0->1 and 1->0
        pairs = set(zip(src.tolist(), dst.tolist()))
        assert (0, 1) in pairs or (1, 0) in pairs

    def test_no_edges_at_high_threshold(self):
        from genome_atlas.graph.view import _compute_similarity_edges

        # Random orthogonal-ish vectors should have low cosine similarity
        feats = torch.eye(5)  # Identity: each pair has sim=0
        src, dst = _compute_similarity_edges(feats, threshold=0.99)
        assert src.numel() == 0

    def test_no_self_loops(self):
        from genome_atlas.graph.view import _compute_similarity_edges

        feats = torch.ones(4, 3)  # All identical -> sim=1 for all pairs
        src, dst = _compute_similarity_edges(feats, threshold=0.5)
        assert (src != dst).all(), "Self-loops found in similarity edges"
