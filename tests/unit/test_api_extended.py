"""Extended unit tests for Atlas API - covers all graph traversal methods.

Covers lines missed by test_api.py:
  - Atlas.__init__ with targets parquet (lines 55-68)
  - Atlas.__init__ with embeddings (line 93)
  - query_protein - graph + protein_map fallback (lines 113-121)
  - systems() with mechanism_bucket filter (lines 123-132)
  - proteins_with_domain (lines 134-147)
  - rna_guides_of_system (lines 162-179)
  - structurally_similar (lines 181-215)
  - structures_of_protein (lines 217-228)
  - get_embedding / similar_nodes (lines 232-257)
  - select_editor (lines 261-274)
  - RuntimeError guards when no graph loaded (lines 126, 138, 172, 192, 220)
"""
from __future__ import annotations

import pickle
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import pytest

from genome_atlas.api import Atlas


# ---------------------------------------------------------------------------
# Helpers - build mock graph and supporting data files
# ---------------------------------------------------------------------------

def _make_full_graph() -> nx.MultiDiGraph:
    G = nx.MultiDiGraph()

    # System nodes
    G.add_node("System_SpCas9", node_type="System", name="SpCas9",
               mechanism_bucket="DSB_NUCLEASE")
    G.add_node("System_PE2",    node_type="System", name="PE2",
               mechanism_bucket="DSB_FREE_PRIME_EDITOR")

    # Protein nodes
    G.add_node("Protein_Q99ZW2", node_type="Protein", accession="Q99ZW2",
               name="Cas9 protein")
    G.add_node("Protein_P12345", node_type="Protein", accession="P12345",
               name="RT subunit")

    # Domain nodes
    G.add_node("Domain_PF09650", node_type="Domain", accession="PF09650",
               name="HNH nuclease domain")
    G.add_node("Domain_PF01535", node_type="Domain", accession="PF01535",
               name="PPR repeat")

    # RNA nodes
    G.add_node("RNA_sgRNA",   node_type="RNA", name="sgRNA",
               rna_type="guide_RNA", length_nt=100)
    G.add_node("RNA_pegRNA",  node_type="RNA", name="pegRNA",
               rna_type="pegRNA",    length_nt=120)

    # Structure nodes
    G.add_node("Structure_7S7W", node_type="Structure", pdb_id="7S7W")
    G.add_node("Structure_6VPC", node_type="Structure", pdb_id="6VPC")

    # Edges - HAS_PROTEIN
    G.add_edge("System_SpCas9", "Protein_Q99ZW2", edge_type="HAS_PROTEIN")
    G.add_edge("System_PE2",    "Protein_P12345", edge_type="HAS_PROTEIN")

    # Edges - HAS_RNA
    G.add_edge("System_SpCas9", "RNA_sgRNA",  edge_type="HAS_RNA")
    G.add_edge("System_PE2",    "RNA_pegRNA", edge_type="HAS_RNA")

    # Edges - HAS_DOMAIN
    G.add_edge("Protein_Q99ZW2", "Domain_PF09650", edge_type="HAS_DOMAIN")
    G.add_edge("Protein_P12345", "Domain_PF01535", edge_type="HAS_DOMAIN")

    # Edges - STRUCTURE_OF (Structure->Protein - used by structures_of_protein)
    G.add_edge("Structure_7S7W", "Protein_Q99ZW2", edge_type="STRUCTURE_OF")

    # Edges - STRUCTURE_OF (Protein->Structure - used by structurally_similar)
    G.add_edge("Protein_Q99ZW2", "Structure_7S7W", edge_type="STRUCTURE_OF")

    # Edges - SIMILAR_TO (Structure->Structure)
    G.add_edge("Structure_7S7W", "Structure_6VPC", edge_type="SIMILAR_TO",
               weight=0.87)

    return G


def _write_graph(path: Path) -> None:
    with open(path, "wb") as f:
        pickle.dump(_make_full_graph(), f)


def _write_targets(path: Path) -> None:
    df = pd.DataFrame([
        {"accession": "Q99ZW2", "length": 1368, "protein_name": "Cas9",
         "organism_name": "Streptococcus pyogenes"},
        {"accession": "P12345", "length":  660, "protein_name": "RT",
         "organism_name": "Moloney murine leukemia virus"},
    ])
    df.to_parquet(path, index=False)


def _write_embeddings(path: Path) -> None:
    # Minimal embeddings parquet: node_id, node_type, embedding columns
    emb = np.random.default_rng(42).random((4, 8)).astype(np.float32)
    rows = [
        {"node_id": "System_SpCas9", "node_type": "System",
         "embedding": emb[0].tolist()},
        {"node_id": "System_PE2",    "node_type": "System",
         "embedding": emb[1].tolist()},
        {"node_id": "Protein_Q99ZW2", "node_type": "Protein",
         "embedding": emb[2].tolist()},
        {"node_id": "Protein_P12345", "node_type": "Protein",
         "embedding": emb[3].tolist()},
    ]
    pd.DataFrame(rows).to_parquet(path, index=False)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def atlas_full(tmp_path_factory):
    d = tmp_path_factory.mktemp("atlas_full")
    gp = d / "atlas.gpickle"
    tp = d / "targets.parquet"
    ep = d / "embeddings.parquet"
    _write_graph(gp)
    _write_targets(tp)
    _write_embeddings(ep)
    return Atlas(gp, ep, tp)


@pytest.fixture(scope="module")
def atlas_graph_only(tmp_path_factory):
    d = tmp_path_factory.mktemp("atlas_graph_only")
    gp = d / "atlas.gpickle"
    _write_graph(gp)
    return Atlas(gp)


@pytest.fixture(scope="module")
def atlas_empty():
    """Atlas with no files - tests RuntimeError guards."""
    return Atlas()


# ---------------------------------------------------------------------------
# Init - targets + embeddings coverage (lines 55-68, 93)
# ---------------------------------------------------------------------------

class TestAtlasInit:

    def test_length_map_populated_from_targets(self, atlas_full):
        assert atlas_full._length_map.get("Protein_Q99ZW2") == 1368
        assert atlas_full._length_map.get("Protein_P12345") == 660

    def test_protein_map_populated(self, atlas_full):
        meta = atlas_full._protein_map.get("Protein_Q99ZW2")
        assert meta is not None
        assert meta["accession"] == "Q99ZW2"
        assert meta["length"] == 1368

    def test_embeddings_loaded(self, atlas_full):
        assert atlas_full._embeddings is not None
        assert "System_SpCas9" in atlas_full._embeddings.index

    def test_system_by_name_cache(self, atlas_full):
        assert "SpCas9" in atlas_full._system_by_name

    def test_protein_by_accession_cache(self, atlas_full):
        assert "Q99ZW2" in atlas_full._protein_by_accession

    def test_rna_by_name_cache(self, atlas_full):
        assert "sgRNA" in atlas_full._rna_by_name

    def test_domain_by_accession_cache(self, atlas_full):
        assert "PF09650" in atlas_full._domain_by_accession


# ---------------------------------------------------------------------------
# query_protein - both graph + protein_map branches (lines 113-121)
# ---------------------------------------------------------------------------

class TestQueryProtein:

    def test_query_protein_in_graph(self, atlas_full):
        result = atlas_full.query_protein("Q99ZW2")
        assert result["node_id"] == "Protein_Q99ZW2"

    def test_query_protein_returns_protein_map_meta(self, atlas_full):
        result = atlas_full.query_protein("Q99ZW2")
        assert result["length"] == 1368

    def test_query_protein_not_in_graph_but_in_map(self, tmp_path):
        """Hit the 'in protein_map but not in graph' branch (lines 119-120)."""
        tp = tmp_path / "targets.parquet"
        _write_targets(tp)
        atlas = Atlas(targets_path=tp)   # no graph
        result = atlas.query_protein("Q99ZW2")
        assert result["accession"] == "Q99ZW2"

    def test_query_protein_missing_raises(self, atlas_full):
        with pytest.raises(KeyError):
            atlas_full.query_protein("ZZZZZZ")


# ---------------------------------------------------------------------------
# systems() (lines 123-132)
# ---------------------------------------------------------------------------

class TestSystems:

    def test_systems_returns_dataframe(self, atlas_full):
        df = atlas_full.systems()
        assert len(df) == 2

    def test_systems_mechanism_bucket_filter(self, atlas_full):
        df = atlas_full.systems(mechanism_bucket="DSB_NUCLEASE")
        assert len(df) == 1
        assert df.iloc[0]["name"] == "SpCas9"

    def test_systems_filter_empty(self, atlas_full):
        df = atlas_full.systems(mechanism_bucket="NONEXISTENT_BUCKET")
        assert len(df) == 0

    def test_systems_no_graph_raises(self, atlas_empty):
        with pytest.raises(RuntimeError, match="not loaded"):
            atlas_empty.systems()


# ---------------------------------------------------------------------------
# proteins_with_domain (lines 134-147)
# ---------------------------------------------------------------------------

class TestProteinsWithDomain:

    def test_proteins_with_domain_found(self, atlas_full):
        df = atlas_full.proteins_with_domain("PF09650")
        assert len(df) == 1
        assert df.iloc[0]["accession"] == "Q99ZW2"

    def test_proteins_with_domain_missing_raises(self, atlas_full):
        with pytest.raises(KeyError):
            atlas_full.proteins_with_domain("PF99999")

    def test_proteins_with_domain_no_graph_raises(self, atlas_empty):
        with pytest.raises(RuntimeError, match="not loaded"):
            atlas_empty.proteins_with_domain("PF09650")


# ---------------------------------------------------------------------------
# rna_guides_of_system (lines 162-179)
# ---------------------------------------------------------------------------

class TestRnaGuidesOfSystem:

    def test_rna_guides_of_spcas9(self, atlas_full):
        df = atlas_full.rna_guides_of_system("SpCas9")
        assert len(df) == 1
        assert df.iloc[0]["name"] == "sgRNA"

    def test_rna_guides_of_pe2(self, atlas_full):
        df = atlas_full.rna_guides_of_system("PE2")
        assert len(df) == 1
        assert df.iloc[0]["name"] == "pegRNA"

    def test_rna_guides_missing_system_raises(self, atlas_full):
        with pytest.raises(KeyError):
            atlas_full.rna_guides_of_system("NotASystem")

    def test_rna_guides_no_graph_raises(self, atlas_empty):
        with pytest.raises(RuntimeError, match="not loaded"):
            atlas_empty.rna_guides_of_system("SpCas9")


# ---------------------------------------------------------------------------
# structures_of_protein (lines 217-228)
# ---------------------------------------------------------------------------

class TestStructuresOfProtein:

    def test_structures_of_q99zw2(self, atlas_full):
        df = atlas_full.structures_of_protein("Q99ZW2")
        assert len(df) == 1
        assert df.iloc[0]["pdb_id"] == "7S7W"

    def test_structures_of_missing_raises(self, atlas_full):
        with pytest.raises(KeyError):
            atlas_full.structures_of_protein("ZZZZZZ")

    def test_structures_no_graph_raises(self, atlas_empty):
        with pytest.raises(RuntimeError, match="not loaded"):
            atlas_empty.structures_of_protein("Q99ZW2")


# ---------------------------------------------------------------------------
# structurally_similar (lines 181-215)
# ---------------------------------------------------------------------------

class TestStructurallySimilar:

    def test_structurally_similar_returns_df(self, atlas_full):
        df = atlas_full.structurally_similar("Q99ZW2", top_k=5)
        # The mock has Protein_Q99ZW2 --STRUCTURE_OF--> Structure_7S7W
        # and Structure_7S7W --SIMILAR_TO--> Structure_6VPC
        assert isinstance(df, pd.DataFrame)

    def test_structurally_similar_tmscore_sorted(self, atlas_full):
        df = atlas_full.structurally_similar("Q99ZW2", top_k=5)
        if len(df) > 1:
            assert list(df["tmscore"]) == sorted(df["tmscore"], reverse=True)

    def test_structurally_similar_missing_raises(self, atlas_full):
        with pytest.raises(KeyError):
            atlas_full.structurally_similar("ZZZZZZ")

    def test_structurally_similar_no_graph_raises(self, atlas_empty):
        with pytest.raises(RuntimeError, match="not loaded"):
            atlas_empty.structurally_similar("Q99ZW2")

    def test_structurally_similar_no_structs_returns_empty(self, atlas_graph_only):
        """P12345 has no Protein->Structure edges so hits the empty-result branch."""
        df = atlas_graph_only.structurally_similar("P12345")
        assert isinstance(df, pd.DataFrame)
        assert df.empty


# ---------------------------------------------------------------------------
# get_embedding / similar_nodes (lines 232-257)
# ---------------------------------------------------------------------------

class TestEmbeddings:

    def test_get_embedding_returns_array(self, atlas_full):
        emb = atlas_full.get_embedding("System_SpCas9")
        assert isinstance(emb, np.ndarray)
        assert emb.shape == (8,)

    def test_get_embedding_missing_raises(self, atlas_full):
        with pytest.raises(KeyError):
            atlas_full.get_embedding("System_NONEXISTENT")

    def test_get_embedding_no_embeddings_raises(self, atlas_graph_only):
        with pytest.raises(RuntimeError, match="No embeddings"):
            atlas_graph_only.get_embedding("System_SpCas9")

    def test_similar_nodes_returns_dataframe(self, atlas_full):
        df = atlas_full.similar_nodes("System_SpCas9", node_type="System", top_k=1)
        assert isinstance(df, pd.DataFrame)

    def test_similar_nodes_empty_type_returns_empty(self, atlas_full):
        df = atlas_full.similar_nodes("System_SpCas9", node_type="NonExistentType")
        assert df.empty

    def test_similar_nodes_no_embeddings_raises(self, atlas_graph_only):
        with pytest.raises(RuntimeError, match="No embeddings"):
            atlas_graph_only.similar_nodes("System_SpCas9")


# ---------------------------------------------------------------------------
# select_editor via Atlas (lines 261-274)
# ---------------------------------------------------------------------------

class TestSelectEditor:

    def test_select_editor_returns_list(self, atlas_full):
        recs = atlas_full.select_editor(
            cell_type="HEK293T", edit_type="deletion",
            cargo_size_bp=0, delivery="AAV", top_k=2
        )
        assert isinstance(recs, list)
        assert len(recs) <= 2

    def test_select_editor_pen_score_range(self, atlas_full):
        recs = atlas_full.select_editor(
            cell_type="HEK293T", edit_type="deletion",
            cargo_size_bp=0, delivery="AAV"
        )
        for r in recs:
            assert 0.0 <= r.pen_score <= 1.0
