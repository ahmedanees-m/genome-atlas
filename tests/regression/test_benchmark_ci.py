"""Regression: benchmark AUROC values must not drift from the v0.6.0 baseline.

JSON key structure (v0.6.0):
    graphsage, gat                  — inductive GNNs (primary benchmark)
    node2vec_inductive              — topology-only, train walks only
    node2vec_transductive           — full-graph walks; supplementary only
    quantum_kernel, classical_rbf   — node classification; supplementary only
"""
import json
from pathlib import Path

import pytest


RESULTS = Path(__file__).parents[2] / "notebooks/benchmark_results.json"


def _load():
    if not RESULTS.exists():
        pytest.skip("benchmark_results.json not found")
    with RESULTS.open() as f:
        return json.load(f)


# ------------------------------------------------------------------ #
# Primary benchmark — inductive GNNs
# ------------------------------------------------------------------ #

def test_graphsage_auroc_above_0_96():
    data = _load()
    auroc = data["graphsage"]["auroc"]
    assert auroc > 0.96, f"GraphSAGE AUROC regressed: {auroc:.4f} <= 0.96"


def test_gat_auroc_above_0_96():
    data = _load()
    auroc = data["gat"]["auroc"]
    assert auroc > 0.96, f"GAT AUROC regressed: {auroc:.4f} <= 0.96"


def test_gat_ci_method_is_bootstrap():
    """v0.6.0 uses 1000× bootstrap CIs, not mann_whitney_se."""
    data = _load()
    method = data["gat"].get("ci_method", "")
    assert "bootstrap" in method, (
        f"GAT CI method should be bootstrap (got '{method}'). "
        "If reverting to mann_whitney_se, update this test and document the reason."
    )


def test_graphsage_ci_method_is_bootstrap():
    data = _load()
    method = data["graphsage"].get("ci_method", "")
    assert "bootstrap" in method, (
        f"GraphSAGE CI method should be bootstrap (got '{method}')."
    )


def test_gat_has_residual_connections():
    """GAT v0.6.0 requires residual connections to prevent embedding collapse."""
    data = _load()
    assert data["gat"].get("residual_connections") is True, (
        "GAT must use residual_connections=True. "
        "Without it, 9,991/10,000 Protein nodes collapse to zero embeddings."
    )


def test_gnns_are_statistically_comparable():
    """GraphSAGE and GAT CIs must overlap — they are tied, not ranked."""
    data = _load()
    sage_lo = data["graphsage"]["ci_lo"]
    sage_hi = data["graphsage"]["ci_hi"]
    gat_lo  = data["gat"]["ci_lo"]
    gat_hi  = data["gat"]["ci_hi"]
    overlap = min(sage_hi, gat_hi) > max(sage_lo, gat_lo)
    assert overlap, (
        f"GNN CIs do not overlap: GraphSAGE [{sage_lo}, {sage_hi}] "
        f"vs GAT [{gat_lo}, {gat_hi}]. "
        "If the gap genuinely widens, remove this test and document the separation."
    )


# ------------------------------------------------------------------ #
# Node2Vec — inductive (topology only)
# ------------------------------------------------------------------ #

def test_node2vec_inductive_auroc_above_0_85():
    """Inductive Node2Vec (train walks only) should comfortably exceed 0.85."""
    data = _load()
    if "node2vec_inductive" not in data:
        pytest.skip("node2vec_inductive not in benchmark_results.json")
    auroc = data["node2vec_inductive"]["auroc"]
    assert auroc > 0.85, (
        f"Inductive Node2Vec AUROC regressed: {auroc:.4f} <= 0.85"
    )


def test_node2vec_excluded_from_primary():
    """Node2Vec entries must be flagged as excluded from the primary table."""
    data = _load()
    for key in ("node2vec_inductive", "node2vec_transductive"):
        if key in data:
            assert data[key].get("excluded_from_primary_table") is True, (
                f"{key} must have excluded_from_primary_table=True. "
                "Node2Vec is not comparable to inductive GNNs."
            )


# ------------------------------------------------------------------ #
# Supplementary methods — sanity bounds only
# ------------------------------------------------------------------ #

def test_quantum_kernel_auroc_above_0_80():
    data = _load()
    auroc = data["quantum_kernel"]["auroc"]
    assert auroc > 0.80, f"Quantum kernel AUROC: {auroc:.4f} <= 0.80"
