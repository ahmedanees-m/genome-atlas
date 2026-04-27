"""Unit tests for selection engine scoring axes."""
import pandas as pd

from genome_atlas.selection import SelectionEngine, UseCaseProfile


class FakeAtlas:
    def __init__(self):
        import networkx as nx
        self._G = nx.MultiDiGraph()
        self._length_map = {}

    def systems(self):
        return pd.DataFrame([
            {"node_id": "System_SpCas9", "name": "SpCas9", "type": "CRISPR-Cas",
             "mechanism_bucket": "DSB_NUCLEASE"},
            {"node_id": "System_Cre", "name": "Cre_recombinase", "type": "Tyrosine_recombinase",
             "mechanism_bucket": "DSB_FREE_TRANSEST_RECOMBINASE"},
            {"node_id": "System_PE2", "name": "PE2", "type": "Prime_editor",
             "mechanism_bucket": "DSB_FREE_TRANSEST_PRIME_EDITOR"},
        ])


def test_rank_orders_by_score():
    engine = SelectionEngine(FakeAtlas())
    recs = engine.rank("HEK293T", "deletion", 0, "AAV", True, top_k=3)
    scores = [r.pen_score for r in recs]
    assert scores == sorted(scores, reverse=True)


def test_nuclease_ranks_high_for_deletion():
    engine = SelectionEngine(FakeAtlas())
    recs = engine.rank("HEK293T", "deletion", 0, "AAV", True, top_k=3)
    # SpCas9 should be top for simple deletion
    assert recs[0].system == "SpCas9"


def test_prime_editor_in_top_3_for_small_insertion():
    engine = SelectionEngine(FakeAtlas())
    recs = engine.rank("HEK293T", "insertion", 100, "AAV", True, top_k=3)
    # PE2 should appear in top-3 for small insertion
    systems = [r.system for r in recs]
    assert "PE2" in systems


def test_dsb_free_ranks_high_for_large_insertion():
    engine = SelectionEngine(FakeAtlas())
    recs = engine.rank("HEK293T", "insertion", 2000, "AAV", True, top_k=3)
    # Cre should be top for large insertion
    assert recs[0].system == "Cre_recombinase"
