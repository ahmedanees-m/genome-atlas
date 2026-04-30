"""Compute 95% bootstrap confidence intervals for v6 GNN benchmark.

Reads the raw per-edge (y_true, y_score) test-set predictions saved by
train_gnn.py (--output-test-preds) and bootstraps AUROC / AUPRC.

For Node2Vec a logistic-regression probe is trained on the train split and
applied to the held-out test split, using the same add_train_val_test_split
seed so the split is identical to the GNN evaluation.

Usage (inside Docker with -v $DATA:/data -v $REPO:/repo -e PYTHONPATH=/repo):
    python3 /repo/scripts/bootstrap_cis_v6.py \
        --graph          /data/graphs/atlas_train.gpickle \
        --node2vec       /data/embeddings/node2vec_v6.parquet \
        --graphsage-preds /data/embeddings/graphsage_v6_test_preds.parquet \
        --gat-preds      /data/embeddings/gat_v6_test_preds.parquet \
        --output         /data/embeddings/bootstrap_cis_v6.parquet
"""
from __future__ import annotations
import argparse
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score

from genome_atlas.graph.build import add_train_val_test_split, build_pyg_hetero


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def bootstrap_ci(y_true: np.ndarray, y_score: np.ndarray, n_boot: int = 1000, seed: int = 42):
    """Return (mean, lo_2.5%, hi_97.5%) AUROC via bootstrap resampling."""
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


def bootstrap_ci_auprc(y_true: np.ndarray, y_score: np.ndarray, n_boot: int = 1000, seed: int = 42):
    rng = np.random.default_rng(seed)
    scores = []
    n = len(y_true)
    for _ in range(n_boot):
        idx = rng.choice(n, n, replace=True)
        if len(np.unique(y_true[idx])) < 2:
            continue
        scores.append(average_precision_score(y_true[idx], y_score[idx]))
    arr = np.array(scores)
    return float(np.mean(arr)), float(np.percentile(arr, 2.5)), float(np.percentile(arr, 97.5))


def node2vec_probe(G, emb_map: dict, edge_type_key: str, seed: int = 42):
    """Train a LR probe on Node2Vec embeddings for one edge type.

    Uses the same deterministic train/test split as train_gnn.py
    (add_train_val_test_split with seed=42, 80/10/10).
    Returns (y_true, y_score) arrays for the test split.

    Negatives are sampled from the *correct node types* (src_type × dst_type)
    to avoid the trivial-classification problem where the LR just learns
    "is the destination a Domain node?" because most random cross-type pairs
    are easy negatives.
    """
    # Parse node types from key: format is always SrcType_REL_DstType where
    # SrcType and DstType are single words (Protein, Domain, Structure, …).
    parts = edge_type_key.split("_")
    src_type = parts[0]   # first token = source node type
    dst_type = parts[-1]  # last token  = destination node type

    # Collect positive edges for this type using graph attributes
    edges = [
        (u, v) for u, v, d in G.edges(data=True)
        if (f"{G.nodes[u].get('node_type','?')}"
            f"_{d.get('edge_type','?')}"
            f"_{G.nodes[v].get('node_type','?')}") == edge_type_key
    ]
    if len(edges) < 20:
        return None, None

    # Split: 80% train / 10% val (unused) / 10% test — matches
    # add_train_val_test_split ordering so train/test sets are comparable
    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(edges))
    n_train = int(len(edges) * 0.8)
    n_val   = int(len(edges) * 0.1)
    train_idx = idx[:n_train]
    test_idx  = idx[n_train + n_val:]

    train_edges = [edges[i] for i in train_idx]
    test_edges  = [edges[i] for i in test_idx]

    # Type-consistent node pools for negative sampling.
    # Negatives are (src_type, dst_type) pairs not already connected,
    # so the LR must learn embedding similarity — not just node type.
    src_nodes = [n for n in G.nodes() if G.nodes[n].get("node_type") == src_type]
    dst_nodes = [n for n in G.nodes() if G.nodes[n].get("node_type") == dst_type]
    if not src_nodes or not dst_nodes:
        return None, None

    dim = len(next(iter(emb_map.values())))

    def build_Xy(pos_edges, n_neg):
        X, y = [], []
        for u, v in pos_edges:
            eu = np.array(emb_map.get(u, [0.0] * dim))
            ev = np.array(emb_map.get(v, [0.0] * dim))
            X.append(np.concatenate([eu, ev]))
            y.append(1)
        neg_rng = np.random.default_rng(seed + 1)
        count, attempts = 0, 0
        while count < n_neg and attempts < n_neg * 200:
            u = src_nodes[neg_rng.integers(len(src_nodes))]
            v = dst_nodes[neg_rng.integers(len(dst_nodes))]
            if not G.has_edge(u, v):
                eu = np.array(emb_map.get(u, [0.0] * dim))
                ev = np.array(emb_map.get(v, [0.0] * dim))
                X.append(np.concatenate([eu, ev]))
                y.append(0)
                count += 1
            attempts += 1
        return np.array(X, dtype=float), np.array(y, dtype=int)

    X_tr, y_tr = build_Xy(train_edges, len(train_edges))
    X_te, y_te = build_Xy(test_edges, len(test_edges))

    if len(np.unique(y_tr)) < 2 or len(np.unique(y_te)) < 2:
        return None, None

    clf = LogisticRegression(max_iter=1000, random_state=seed)
    clf.fit(X_tr, y_tr)
    y_score = clf.predict_proba(X_te)[:, 1]
    return y_te, y_score


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--graph",            type=Path, required=True)
    p.add_argument("--node2vec",         type=Path, required=True)
    p.add_argument("--graphsage-preds",  type=Path, required=True)
    p.add_argument("--gat-preds",        type=Path, required=True)
    p.add_argument("--output",           type=Path, required=True)
    p.add_argument("--n-boot",           type=int, default=1000)
    args = p.parse_args()

    print("Loading graph...")
    with open(args.graph, "rb") as f:
        G = pickle.load(f)

    print("Loading Node2Vec embeddings...")
    n2v_df = pd.read_parquet(args.node2vec)
    emb_map = dict(zip(n2v_df["node_id"], n2v_df["embedding"].apply(np.array)))

    results = []

    # ------------------------------------------------------------------ #
    # GNN methods: use saved test predictions                             #
    # ------------------------------------------------------------------ #
    for method, preds_path in [("GraphSAGE", args.graphsage_preds), ("GAT", args.gat_preds)]:
        print(f"\n  {method}...")
        df = pd.read_parquet(preds_path)
        for et, grp in df.groupby("edge_type"):
            y_true  = grp["y_true"].values
            y_score = grp["y_score"].values
            if len(np.unique(y_true)) < 2 or len(y_true) < 10:
                print(f"    {et}: skipped (n={len(y_true)}, classes={np.unique(y_true)})")
                continue
            auroc_m, auroc_lo, auroc_hi = bootstrap_ci(y_true, y_score, args.n_boot)
            auprc_m, auprc_lo, auprc_hi = bootstrap_ci_auprc(y_true, y_score, args.n_boot)
            print(f"    {et}: AUROC={auroc_m:.4f} [{auroc_lo:.4f}, {auroc_hi:.4f}]")
            results.append({
                "method": method,
                "edge_type": et,
                "auroc_mean": auroc_m, "auroc_lo": auroc_lo, "auroc_hi": auroc_hi,
                "auprc_mean": auprc_m, "auprc_lo": auprc_lo, "auprc_hi": auprc_hi,
                "n_test": len(y_true),
                # Evaluation mode metadata — critical for manuscript reporting
                "evaluation_mode": "Inductive",
                "notes": "Test edges withheld from message-passing during GNN training",
            })

    # ------------------------------------------------------------------ #
    # Node2Vec: probe per edge type                                       #
    # ------------------------------------------------------------------ #
    print("\n  Node2Vec...")
    edge_types_in_graph = set()
    for u, v, d in G.edges(data=True):
        et = f"{G.nodes[u].get('node_type','?')}_{d.get('edge_type','?')}_{G.nodes[v].get('node_type','?')}"
        edge_types_in_graph.add(et)

    for et in sorted(edge_types_in_graph):
        y_true, y_score = node2vec_probe(G, emb_map, et)
        if y_true is None:
            continue
        if len(np.unique(y_true)) < 2:
            continue
        auroc_m, auroc_lo, auroc_hi = bootstrap_ci(y_true, y_score, args.n_boot)
        auprc_m, auprc_lo, auprc_hi = bootstrap_ci_auprc(y_true, y_score, args.n_boot)
        print(f"    {et}: AUROC={auroc_m:.4f} [{auroc_lo:.4f}, {auroc_hi:.4f}]")
        results.append({
            "method": "Node2Vec",
            "edge_type": et,
            "auroc_mean": auroc_m, "auroc_lo": auroc_lo, "auroc_hi": auroc_hi,
            "auprc_mean": auprc_m, "auprc_lo": auprc_lo, "auprc_hi": auprc_hi,
            "n_test": len(y_true),
            # IMPORTANT: Node2Vec is TRANSDUCTIVE — random walks are trained on the
            # full graph including test edges.  This inflates AUROC and makes it
            # incomparable to inductive GNN results.  Do NOT include Node2Vec in
            # the same primary-benchmark table as GraphSAGE/GAT.
            # Report separately as a topology upper-bound in supplementary.
            "evaluation_mode": "Transductive",
            "notes": "Random walks trained on full graph including test edges; "
                     "upper bound only — not comparable to inductive GNNs",
        })

    # ------------------------------------------------------------------ #
    # Save                                                                #
    # ------------------------------------------------------------------ #
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out_df = pd.DataFrame(results)
    out_df.to_parquet(args.output, compression="zstd")
    print(f"\nSaved bootstrap CIs -> {args.output}")

    # Pretty-print primary benchmark
    print("\n" + "=" * 80)
    print("BOOTSTRAP 95% CIs — PRIMARY BENCHMARK (Protein_HAS_DOMAIN_Domain)")
    print("=" * 80)
    primary = out_df[out_df["edge_type"] == "Protein_HAS_DOMAIN_Domain"]
    for _, row in primary.iterrows():
        print(
            f"  {row['method']:12s}: AUROC = {row['auroc_mean']:.4f} "
            f"[{row['auroc_lo']:.4f}, {row['auroc_hi']:.4f}]  "
            f"AUPRC = {row['auprc_mean']:.4f} "
            f"[{row['auprc_lo']:.4f}, {row['auprc_hi']:.4f}]  "
            f"n={row['n_test']}"
        )
    print("=" * 80)
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
