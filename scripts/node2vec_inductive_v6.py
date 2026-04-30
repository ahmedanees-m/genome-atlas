"""Inductive Node2Vec evaluation for the v6 benchmark.

Trains Node2Vec ONLY on train edges (test and val edges withheld from
random walks), then evaluates link prediction on the held-out test split
using a logistic-regression probe.  This gives a fair, apples-to-apples
comparison with the inductive GNNs.

The exact same 80/10/10 split (seed=42) as train_gnn.py is used by
delegating to add_train_val_test_split, so test-edge sets are identical.

Negatives are sampled from the correct node-type pools (type-consistent),
matching the approach in bootstrap_cis_v6.py.

Output parquet columns:
    method, edge_type, auroc_mean, auroc_lo, auroc_hi,
    auprc_mean, n_test, evaluation_mode, notes

Usage (inside Docker with -v $DATA:/data -v $REPO:/repo -e PYTHONPATH=/repo):
    python3 /repo/scripts/node2vec_inductive_v6.py \\
        --graph          /data/graphs/atlas_train.gpickle \\
        --esm-embeddings /data/embeddings/esm2_150M_v6.parquet \\
        --output         /data/embeddings/node2vec_inductive_v6.parquet
"""
from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score

from genome_atlas.graph.build import add_train_val_test_split, build_pyg_hetero


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def bootstrap_ci(y_true: np.ndarray, y_score: np.ndarray,
                 n_boot: int = 1000, seed: int = 42):
    rng = np.random.default_rng(seed)
    scores = []
    n = len(y_true)
    for _ in range(n_boot):
        idx = rng.choice(n, n, replace=True)
        if len(np.unique(y_true[idx])) < 2:
            continue
        scores.append(roc_auc_score(y_true[idx], y_score[idx]))
    arr = np.array(scores)
    return float(np.mean(arr)), float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))


def build_Xy(pos_edges, src_nodes, dst_nodes, emb_map, dim, rng_seed, G):
    """Build (X, y) from positive edges + type-consistent negatives.

    Args:
        pos_edges:  List of (u, v) positive pairs.
        src_nodes:  All nodes of the source type (for negative sampling).
        dst_nodes:  All nodes of the destination type.
        emb_map:    {node_id: embedding_array}.
        dim:        Embedding dimension (used for missing-node fallback).
        rng_seed:   RNG seed for reproducible negative sampling.
        G:          Full NetworkX graph (used to exclude existing edges from negatives).
    """
    X, y = [], []
    for u, v in pos_edges:
        eu = np.array(emb_map.get(u, [0.0] * dim))
        ev = np.array(emb_map.get(v, [0.0] * dim))
        X.append(np.concatenate([eu, ev]))
        y.append(1)
    neg_rng = np.random.default_rng(rng_seed)
    count = attempts = 0
    n_neg = len(pos_edges)
    while count < n_neg and attempts < n_neg * 300:
        u = src_nodes[neg_rng.integers(len(src_nodes))]
        v = dst_nodes[neg_rng.integers(len(dst_nodes))]
        if not G.has_edge(u, v):
            eu = np.array(emb_map.get(u, [0.0] * dim))
            ev = np.array(emb_map.get(v, [0.0] * dim))
            X.append(np.concatenate([eu, ev]))
            y.append(0)
            count += 1
        attempts += 1
    return np.array(X, dtype=np.float64), np.array(y, dtype=np.int32)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--graph",          type=Path, required=True)
    p.add_argument("--esm-embeddings", type=Path, required=True)
    p.add_argument("--output",         type=Path, required=True)
    p.add_argument("--dimensions",     type=int,  default=128)
    p.add_argument("--walk-length",    type=int,  default=30)
    p.add_argument("--num-walks",      type=int,  default=10)
    p.add_argument("--workers",        type=int,  default=8)
    p.add_argument("--n-boot",         type=int,  default=1000)
    args = p.parse_args()

    # ------------------------------------------------------------------ #
    # Load full graph (needed for negative-edge checks)
    # ------------------------------------------------------------------ #
    print("Loading graph...")
    with open(args.graph, "rb") as f:
        G = pickle.load(f)
    print(f"  {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")

    # ------------------------------------------------------------------ #
    # Build PyG graph → get EXACT same train/test split as train_gnn.py
    # ------------------------------------------------------------------ #
    print("Building PyG graph and applying 80/10/10 split (seed=42)...")
    data = build_pyg_hetero(args.graph, args.esm_embeddings)
    data = add_train_val_test_split(data, train_ratio=0.8, val_ratio=0.1, seed=42)

    # Convert PyG test+val indices back to NetworkX node-ID pairs so we can
    # remove them from the graph before training Node2Vec
    withheld: set[tuple] = set()
    for edge_type in data.edge_types:
        src_type, rel, dst_type = edge_type
        src_ids = data[src_type].node_ids
        dst_ids = data[dst_type].node_ids
        for mask_attr in ("test_mask", "val_mask"):
            mask = getattr(data[edge_type], mask_attr)
            ei   = data[edge_type].edge_index[:, mask]
            for s, d in zip(ei[0].tolist(), ei[1].tolist()):
                withheld.add((src_ids[s], dst_ids[d]))

    print(f"  Withholding {len(withheld):,} test+val edges from Node2Vec walks")

    # ------------------------------------------------------------------ #
    # Build train-only NetworkX graph
    # ------------------------------------------------------------------ #
    G_train = nx.DiGraph()
    G_train.add_nodes_from(G.nodes(data=True))
    kept = skipped = 0
    for u, v, d in G.edges(data=True):
        if (u, v) in withheld:
            skipped += 1
        else:
            G_train.add_edge(u, v, **d)
            kept += 1
    print(f"  Train graph: {G_train.number_of_nodes():,} nodes, "
          f"{G_train.number_of_edges():,} edges  (removed {skipped:,})")

    # ------------------------------------------------------------------ #
    # Train Node2Vec on train-only graph
    # ------------------------------------------------------------------ #
    print(f"\nTraining Node2Vec  dim={args.dimensions}  "
          f"walks={args.num_walks}×{args.walk_length}  workers={args.workers}...")
    from node2vec import Node2Vec as _Node2Vec
    n2v   = _Node2Vec(
        G_train,
        dimensions=args.dimensions,
        walk_length=args.walk_length,
        num_walks=args.num_walks,
        workers=args.workers,
        seed=42,
        quiet=False,
    )
    model   = n2v.fit(window=10, min_count=1, batch_words=4)
    emb_map = {node: model.wv[node] for node in G_train.nodes() if node in model.wv}
    print(f"  Embeddings for {len(emb_map):,} nodes")

    # ------------------------------------------------------------------ #
    # Evaluate per edge type using the held-out TEST split
    # ------------------------------------------------------------------ #
    results = []
    dim     = args.dimensions

    for edge_type in data.edge_types:
        src_type, rel, dst_type = edge_type
        et_str = f"{src_type}_{rel}_{dst_type}"

        src_ids = data[src_type].node_ids
        dst_ids = data[dst_type].node_ids

        test_ei  = data[edge_type].edge_index[:, data[edge_type].test_mask]
        train_ei = data[edge_type].edge_index[:, data[edge_type].train_mask]

        if test_ei.size(1) < 10:
            print(f"  {et_str}: skipped (n_test={test_ei.size(1)})")
            continue

        test_pos  = [(src_ids[s], dst_ids[d])
                     for s, d in zip(test_ei[0].tolist(),  test_ei[1].tolist())]
        train_pos = [(src_ids[s], dst_ids[d])
                     for s, d in zip(train_ei[0].tolist(), train_ei[1].tolist())]

        # Type-consistent node pools for negative sampling
        src_nodes = [n for n in G.nodes()
                     if G.nodes[n].get("node_type") == src_type]
        dst_nodes = [n for n in G.nodes()
                     if G.nodes[n].get("node_type") == dst_type]
        if not src_nodes or not dst_nodes:
            continue

        X_tr, y_tr = build_Xy(train_pos, src_nodes, dst_nodes, emb_map, dim, rng_seed=42, G=G)
        X_te, y_te = build_Xy(test_pos,  src_nodes, dst_nodes, emb_map, dim, rng_seed=43, G=G)

        if len(np.unique(y_tr)) < 2 or len(np.unique(y_te)) < 2:
            print(f"  {et_str}: skipped (single class in train or test)")
            continue

        clf = LogisticRegression(max_iter=1000, random_state=42)
        clf.fit(X_tr, y_tr)
        y_score = clf.predict_proba(X_te)[:, 1]

        auroc_m, auroc_lo, auroc_hi = bootstrap_ci(y_te, y_score, args.n_boot)
        auprc_m = float(average_precision_score(y_te, y_score))

        print(f"  {et_str}: AUROC={auroc_m:.4f} [{auroc_lo:.4f}, {auroc_hi:.4f}]"
              f"  AUPRC={auprc_m:.4f}  n={len(y_te)}")

        results.append({
            "method":          "Node2Vec",
            "edge_type":       et_str,
            "auroc_mean":      auroc_m,
            "auroc_lo":        auroc_lo,
            "auroc_hi":        auroc_hi,
            "auprc_mean":      auprc_m,
            "auprc_lo":        None,
            "auprc_hi":        None,
            "n_test":          len(y_te),
            "evaluation_mode": "Inductive",
            "notes":           "Node2Vec trained on train edges only; "
                               "test/val edges withheld from random walks",
        })

    # ------------------------------------------------------------------ #
    # Save
    # ------------------------------------------------------------------ #
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(results)
    df.to_parquet(args.output, compression="zstd")
    print(f"\nSaved -> {args.output}  ({len(df)} rows)")

    # Print primary benchmark summary
    primary = df[df["edge_type"] == "Protein_HAS_DOMAIN_Domain"]
    if len(primary):
        row = primary.iloc[0]
        print("\n" + "=" * 72)
        print("INDUCTIVE Node2Vec — Protein_HAS_DOMAIN_Domain")
        print("=" * 72)
        print(f"  AUROC = {row['auroc_mean']:.4f}  [{row['auroc_lo']:.4f}, {row['auroc_hi']:.4f}]")
        print(f"  AUPRC = {row['auprc_mean']:.4f}")
        print(f"  n_test = {row['n_test']}")
        print(f"  Transductive AUROC = 0.9965  →  Δ = {0.9965 - row['auroc_mean']:.4f}")
        print("=" * 72)

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
