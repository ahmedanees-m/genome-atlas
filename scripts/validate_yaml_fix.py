"""Quick validation of engineered variants without full graph."""
import yaml
from pathlib import Path

import pandas as pd

from genome_atlas.selection import SelectionEngine


def main():
    # Load scenarios
    scenarios_path = Path("notebooks/validation_scenarios.yaml")
    with scenarios_path.open() as f:
        data = yaml.safe_load(f)

    # Build mock atlas from YAML systems
    systems_path = Path("genome_atlas/data/foundational_systems.yaml")
    with systems_path.open() as fy:
        sys_data = yaml.safe_load(fy)

    import networkx as nx

    class MockAtlas:
        def __init__(self):
            self._G = nx.MultiDiGraph()
            self._length_map = {}

        def systems(self):
            rows = []
            for s in sys_data["systems"]:
                rows.append({
                    "node_id": f"System_{s['name']}",
                    "name": s["name"],
                    "type": s["type"],
                    "mechanism_bucket": s["mechanism_bucket"],
                })
            return pd.DataFrame(rows)

    atlas = MockAtlas()
    engine = SelectionEngine(atlas)

    hits = 0
    total = 0
    for s in data["scenarios"]:
        p = s["params"]
        recs = engine.rank(
            p["cell_type"], p["edit_type"], p["cargo_size_bp"],
            p["delivery"], True, top_k=5,
        )
        top3 = [r.system for r in recs[:3]]
        published = s["published_editor"]
        ok = published in top3
        hits += ok
        total += 1
        flag = "PASS" if ok else "FAIL"
        print(f"{flag} {s['name']}: published={published} | top3={top3}")

    print(f"\nAccuracy: {hits}/{total} = {100*hits/total:.1f}%")


if __name__ == "__main__":
    main()
