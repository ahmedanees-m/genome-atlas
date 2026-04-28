"""Diagnostic: check system_total_size_aa for all systems."""
import pickle
from pathlib import Path

import networkx as nx
import yaml

from genome_atlas.utils.size import system_total_size_aa


def main():
    graph_path = Path.home() / "pen-stack/data/graphs/atlas.gpickle"
    if graph_path.exists():
        with open(graph_path, "rb") as f:
            G = pickle.load(f)
        print("Loaded real graph")
    else:
        G = nx.MultiDiGraph()
        print("Real graph not found — using empty mock")

    class FakeAtlas:
        def __init__(self):
            self._G = G
            self._length_map = {}

    atlas = FakeAtlas()

    with open("genome_atlas/data/foundational_systems.yaml") as f:
        data = yaml.safe_load(f)

    print(f"{'System':30s} {'total_aa':>8s}  -> AAV score")
    print("-" * 55)
    for s in data["systems"]:
        nid = f"System_{s['name']}"
        total = system_total_size_aa(atlas, nid)
        if total == 0:
            aav = 0.1  # unknown size guard
        else:
            aav = 1.0 if total <= 900 else (0.6 if total <= 2000 else 0.2)
        print(f"{s['name']:30s} {total:8d} aa  -> {aav:.1f}")


if __name__ == "__main__":
    main()
