"""Regression: benchmark AUROC values must not drift from recorded baseline."""
import json
from pathlib import Path

import pytest


RESULTS = Path(__file__).parents[2] / "notebooks/benchmark_results.json"


def _load():
    if not RESULTS.exists():
        pytest.skip("benchmark_results.json not found")
    with RESULTS.open() as f:
        return json.load(f)


def test_gat_auroc_above_0_96():
    data = _load()
    assert data["gat"]["auroc"] > 0.96, (
        f"GAT AUROC regressed: {data['gat']['auroc']:.4f} <= 0.96"
    )


def test_graphsage_auroc_above_0_95():
    data = _load()
    assert data["graphsage"]["auroc"] > 0.95, (
        f"GraphSAGE AUROC regressed: {data['graphsage']['auroc']:.4f} <= 0.95"
    )


def test_node2vec_auroc_above_0_80():
    data = _load()
    assert data["node2vec"]["auroc"] > 0.80, (
        f"Node2Vec AUROC regressed: {data['node2vec']['auroc']:.4f} <= 0.80"
    )


def test_gat_ci_method_is_mann_whitney():
    data = _load()
    assert data["gat"].get("ci_method") == "mann_whitney_se", (
        "GAT CI method must be mann_whitney_se to distinguish from inflated bootstrap"
    )
