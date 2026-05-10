"""Shared test fixtures."""
from pathlib import Path
import pytest


# VM data paths (used when running on the compute server)
_ATLAS_BASE = Path.home() / "pen-stack/data"


@pytest.fixture
def demo_graph_path():
    p = _ATLAS_BASE / "graphs/atlas.gpickle"
    if not p.exists():
        pytest.skip("Full atlas not available on this runner")
    return p


@pytest.fixture
def demo_embeddings_path():
    return _ATLAS_BASE / "embeddings/graphsage.parquet"


@pytest.fixture
def demo_targets_path():
    p = _ATLAS_BASE / "processed/targets_v2.parquet"
    if not p.exists():
        pytest.skip("targets_v2 not available on this runner")
    return p
