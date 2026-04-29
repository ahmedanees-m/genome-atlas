"""Create atlas_train.gpickle: full atlas with isolated nodes removed.

Isolated nodes (degree = 0) are Structure or Domain orphans that have neither
a STRUCTURE_OF edge to any protein in the atlas nor a SIMILAR_TO edge to any
other structure.  They contribute nothing to GNN message passing or Node2Vec
random walks, but inflate node counts and dilute embedding quality.

The full atlas.gpickle is kept intact for the API. This training variant is
used exclusively by the GNN / Node2Vec training pipeline.

Usage (on VM inside pen-stack/graph:0.1.0 Docker):
    python3 scripts/create_train_gpickle.py
    # or with explicit paths:
    python3 scripts/create_train_gpickle.py \
        --input  /data/graphs/atlas.gpickle \
        --output /data/graphs/atlas_train.gpickle
"""
import argparse
import os
import pickle
from pathlib import Path

import networkx as nx


def main(input_path: Path, output_path: Path) -> None:
    print("=" * 60)
    print("Create Training GPicle (remove isolated nodes)")
    print("=" * 60)

    print(f"  Loading: {input_path}")
    with open(str(input_path), "rb") as f:
        G = pickle.load(f)

    n_before = G.number_of_nodes()
    e_before = G.number_of_edges()
    print(f"  Before: {n_before:,} nodes, {e_before:,} edges")

    isolates = list(nx.isolates(G))
    type_counts: dict[str, int] = {}
    for nid in isolates:
        nt = G.nodes[nid].get("node_type", "Unknown")
        type_counts[nt] = type_counts.get(nt, 0) + 1
    print(f"  Isolated nodes to remove: {len(isolates)}")
    for nt, cnt in sorted(type_counts.items(), key=lambda x: -x[1]):
        print(f"    {nt}: {cnt}")

    G.remove_nodes_from(isolates)

    n_after = G.number_of_nodes()
    e_after = G.number_of_edges()
    print(f"  After:  {n_after:,} nodes, {e_after:,} edges")

    # Verify no remaining isolates
    remaining = len(list(nx.isolates(G)))
    print(f"  Remaining isolates: {remaining} (should be 0)")

    # Component stats
    comps = sorted(
        [len(c) for c in nx.weakly_connected_components(G)], reverse=True
    )
    print(f"  Connected components: {len(comps)}")
    print(f"  Largest 5: {comps[:5]}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(str(output_path), "wb") as f:
        pickle.dump(G, f)
    sz = output_path.stat().st_size / (1024 * 1024)
    print(f"  Saved: {output_path}  ({sz:.1f} MB)")
    print("\nDone. Next: bash pipeline_v6.sh")


if __name__ == "__main__":
    _base = Path(os.environ.get("PENSTACK_DATA", "/data"))
    ap = argparse.ArgumentParser()
    ap.add_argument("--input",  type=Path, default=_base / "graphs/atlas.gpickle")
    ap.add_argument("--output", type=Path, default=_base / "graphs/atlas_train.gpickle")
    args = ap.parse_args()
    main(args.input, args.output)
