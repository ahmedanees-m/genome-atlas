"""Validation of selection decision support against published use cases."""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).parent.parent))
from genome_atlas.api import Atlas


def main():
    atlas = Atlas.load(
        gpickle_path=Path("/data/graphs/atlas.gpickle"),
        embeddings_path=Path("/data/embeddings/graphsage.parquet"),
        targets_path=Path("/data/processed/targets_v2.parquet"),
    )

    cfg_path = Path(__file__).with_suffix("").with_suffix(".yaml")
    if not cfg_path.exists():
        cfg_path = Path(__file__).parent / "validation_scenarios.yaml"

    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    rows = []
    for s in cfg["scenarios"]:
        recs = atlas.select_editor(**s["params"], top_k=5)
        top_names = [r.system for r in recs]
        hit = s["published_editor"] in top_names[:3]
        rank = (top_names.index(s["published_editor"]) + 1
                if s["published_editor"] in top_names else None)
        rows.append({
            "scenario": s["name"],
            "published_editor": s["published_editor"],
            "top3_recommendations": top_names[:3],
            "published_in_top3": hit,
            "published_rank": rank,
        })

    df = pd.DataFrame(rows)
    print(df.to_string(index=False))
    accuracy = df["published_in_top3"].mean()
    print(f"\nTop-3 accuracy: {accuracy:.1%}")

    # Save results
    out_path = Path("/data/embeddings/selection_validation.parquet")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out_path, compression="zstd")
    print(f"Saved -> {out_path}")

    assert accuracy >= 0.7, f"Validation gate failed: top-3 accuracy {accuracy:.1%} < 70%"
    print("Validation gate PASSED.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
