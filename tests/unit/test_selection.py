"""Unit tests for selection engine."""
import pytest

from genome_atlas.selection import SelectionEngine, UseCaseProfile


class FakeAtlas:
    def __init__(self):
        import networkx as nx
        self._G = nx.MultiDiGraph()
        self._length_map = {}

    def systems(self):
        import pandas as pd
        return pd.DataFrame([
            {"node_id": "System_SpCas9", "name": "SpCas9", "type": "CRISPR-Cas",
             "mechanism_bucket": "DSB_NUCLEASE"},
            {"node_id": "System_Cre", "name": "Cre_recombinase", "type": "Tyrosine_recombinase",
             "mechanism_bucket": "DSB_FREE_TRANSEST_RECOMBINASE"},
        ])


def test_rank_returns_top_k():
    engine = SelectionEngine(FakeAtlas())
    recs = engine.rank("HEK293T", "deletion", 0, "AAV", True, top_k=2)
    assert len(recs) == 2
    assert recs[0].system == "SpCas9" or recs[0].system == "Cre_recombinase"


def test_no_match_guard():
    class EmptyAtlas:
        def systems(self):
            import pandas as pd
            return pd.DataFrame()

    engine = SelectionEngine(EmptyAtlas())
    recs = engine.rank("HEK293T", "deletion", 0, "AAV", True)
    assert len(recs) == 1
    assert recs[0].system == "NO_MATCH"
