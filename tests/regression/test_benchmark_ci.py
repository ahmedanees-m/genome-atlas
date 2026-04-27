"""Regression: benchmark values must not drift."""
import json
from pathlib import Path

import pytest


RESULTS = Path(__file__).parents[3] / "notebooks/benchmark_results.json"


def _load():
    if not RESULTS.exists():
        pytest.skip("benchmark_results.json not found")
    with RESULTS.open() as f:
        return json.load(f)


def test_gat_auroc_above_0_96():
    data = _load()
    assert data["gat"]["auroc"] > 0.96


def test_graphsage_auroc_above_0_95():
    data = _load()
    assert data["graphsage"]["auroc"] > 0.95


def test_node2vec_auroc_above_0_80():
    data = _load()
    assert data["node2vec"]["auroc"] > 0.80
