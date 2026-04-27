"""Compute total coding-sequence size of a system (sum of its protein subunits)."""
from __future__ import annotations

from genome_atlas.api import Atlas


def system_total_size_aa(atlas: Atlas, system_node_id: str) -> int:
    """Sum protein lengths for all proteins linked to a system via HAS_PROTEIN."""
    if atlas._G is None:
        return 0
    total = 0
    for u, v, d in atlas._G.edges(data=True):
        if d.get("edge_type") == "HAS_PROTEIN" and u == system_node_id:
            total += atlas._length_map.get(v, 0)
    return total
