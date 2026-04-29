"""Train heterogeneous GNN (GraphSAGE or GAT) with link prediction.

Differences from archive/train_gnn.py:
- in_features is inferred from data.x_dict (no hardcoded 480 / 640).
- evaluate() uses a fixed RNG seed for reproducible negative sampling.
- Saves raw test-set (y_true, y_score) pairs to --output-test-preds for
  downstream bootstrap CI computation.
"""
from __future__ import annotations
import argparse
import pickle
from pathlib import Path
from typing import Tuple

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from sklearn.metrics import roc_auc_score, average_precision_score

from genome_atlas.graph.build import build_pyg_hetero, add_train_val_test_split
from genome_atlas.models.graphsage import HeteroGNN, LinkPredictor


def compute_loss(model, predictor, data, edge_type, mask):
    src_type, rel, dst_type = edge_type
    edge_index = data[edge_type].edge_index[:, mask]
    if edge_index.size(1) == 0:
        return torch.tensor(0.0, device=next(model.parameters()).device)

    out = model(data.x_dict, data.edge_index_dict)
    src_emb = out[src_type][edge_index[0]]
    dst_emb = out[dst_type][edge_index[1]]
    pos_score = predictor(src_emb, dst_emb)

    num_neg = edge_index.size(1)
    num_dst = data[dst_type].x.size(0)
    neg_dst = torch.randint(0, num_dst, (num_neg,), device=edge_index.device)
    neg_emb = out[dst_type][neg_dst]
    neg_score = predictor(src_emb, neg_emb)

    scores = torch.cat([pos_score, neg_score])
    labels = torch.cat([torch.ones_like(pos_score), torch.zeros_like(neg_score)])
    return F.binary_cross_entropy(scores, labels)


def evaluate(
    model, predictor, data, edge_type, mask, seed: int = 42
) -> Tuple[dict, np.ndarray, np.ndarray]:
    """Return (metrics_dict, y_true, y_score) with fixed-seed negatives."""
    model.eval()
    src_type, rel, dst_type = edge_type
    edge_index = data[edge_type].edge_index[:, mask]
    if edge_index.size(1) == 0:
        return {"auroc": 0.0, "auprc": 0.0}, np.array([]), np.array([])

    with torch.no_grad():
        out = model(data.x_dict, data.edge_index_dict)
        src_emb = out[src_type][edge_index[0]]
        dst_emb = out[dst_type][edge_index[1]]
        pos_score = predictor(src_emb, dst_emb).cpu().numpy()

        num_neg = edge_index.size(1)
        num_dst = data[dst_type].x.size(0)
        # Fixed seed: reproducible negatives so test AUROC is deterministic
        gen = torch.Generator(device=edge_index.device)
        gen.manual_seed(seed)
        neg_dst = torch.randint(0, num_dst, (num_neg,), device=edge_index.device, generator=gen)
        neg_emb = out[dst_type][neg_dst]
        neg_score = predictor(src_emb, neg_emb).cpu().numpy()

    y_score = np.concatenate([pos_score, neg_score])
    y_true  = np.concatenate([np.ones_like(pos_score), np.zeros_like(neg_score)])
    auroc = roc_auc_score(y_true, y_score)
    auprc = average_precision_score(y_true, y_score)
    return {"auroc": auroc, "auprc": auprc}, y_true, y_score


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--model", choices=["sage", "gat"], required=True)
    p.add_argument("--graph", type=Path, required=True)
    p.add_argument("--esm-embeddings", type=Path, required=True)
    p.add_argument("--output-embeddings", type=Path, required=True)
    p.add_argument("--output-metrics", type=Path, required=True)
    p.add_argument("--output-test-preds", type=Path, default=None,
                   help="Parquet with (edge_type, y_true, y_score) for bootstrap CIs")
    p.add_argument("--epochs", type=int, default=100)
    p.add_argument("--lr", type=float, default=0.001)
    p.add_argument("--hidden", type=int, default=128)
    p.add_argument("--out-dim", type=int, default=128)
    p.add_argument("--weight-decay", type=float, default=1e-4)
    args = p.parse_args()

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    print("Building PyG graph...")
    data = build_pyg_hetero(args.graph, args.esm_embeddings)
    data = add_train_val_test_split(data)
    data = data.to(device)
    print(f"Node types: {list(data.node_types)}")
    print(f"Edge types: {list(data.edge_types)}")

    # Infer input feature dimension from data (handles any ESM-2 variant)
    in_features = next(iter(data.x_dict.values())).size(1)
    print(f"Inferred in_features={in_features} from x_dict")

    model = HeteroGNN(
        data.metadata(),
        hidden_channels=args.hidden,
        out_channels=args.out_dim,
        model_type=args.model,
        in_features=in_features,
    ).to(device)
    predictor = LinkPredictor(args.out_dim).to(device)
    optimizer = torch.optim.Adam(
        list(model.parameters()) + list(predictor.parameters()),
        lr=args.lr, weight_decay=args.weight_decay,
    )

    best_val_auroc = 0.0
    best_state = None

    print(f"Training {args.model.upper()} for {args.epochs} epochs...")
    for epoch in range(1, args.epochs + 1):
        model.train()
        total_loss = 0.0
        for edge_type in data.edge_types:
            mask = data[edge_type].train_mask
            if mask.sum() == 0:
                continue
            optimizer.zero_grad()
            loss = compute_loss(model, predictor, data, edge_type, mask)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()

        if epoch % 10 == 0 or epoch == 1:
            val_auroc_sum = 0.0
            val_count = 0
            for edge_type in data.edge_types:
                mask = data[edge_type].val_mask
                if mask.sum() == 0:
                    continue
                res, _, _ = evaluate(model, predictor, data, edge_type, mask)
                val_auroc_sum += res["auroc"]
                val_count += 1
            val_auroc = val_auroc_sum / max(1, val_count)
            print(f"  Epoch {epoch:>3d}: loss={total_loss:.4f}, val_auroc={val_auroc:.4f}")
            if val_auroc > best_val_auroc:
                best_val_auroc = val_auroc
                best_state = {
                    "model": {k: v.cpu() for k, v in model.state_dict().items()},
                    "predictor": {k: v.cpu() for k, v in predictor.state_dict().items()},
                }

    # Load best checkpoint
    if best_state:
        model.load_state_dict({k: v.to(device) for k, v in best_state["model"].items()})
        predictor.load_state_dict({k: v.to(device) for k, v in best_state["predictor"].items()})

    # Final evaluation on test set
    print("\nFinal test evaluation:")
    metrics = []
    test_preds_records = []

    for edge_type in data.edge_types:
        src, rel, dst = edge_type
        mask = data[edge_type].test_mask
        if mask.sum() == 0:
            continue
        res, y_true, y_score = evaluate(model, predictor, data, edge_type, mask, seed=42)
        print(f"  {src}->{rel}->{dst}: AUROC={res['auroc']:.4f}, AUPRC={res['auprc']:.4f}")
        et_str = f"{src}_{rel}_{dst}"
        metrics.append({
            "edge_type": et_str,
            "model": args.model,
            "auroc": res["auroc"],
            "auprc": res["auprc"],
        })
        for yt, ys in zip(y_true.tolist(), y_score.tolist()):
            test_preds_records.append({
                "edge_type": et_str,
                "y_true": float(yt),
                "y_score": float(ys),
            })

    # Save node embeddings
    model.eval()
    with torch.no_grad():
        out = model(data.x_dict, data.edge_index_dict)
    records = []
    for node_type in data.node_types:
        for i, node_id in enumerate(data[node_type].node_ids):
            records.append({
                "node_id": node_id,
                "node_type": node_type,
                "embedding": out[node_type][i].cpu().numpy().tolist(),
            })
    args.output_embeddings.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(records).to_parquet(args.output_embeddings, compression="zstd")
    print(f"Saved embeddings -> {args.output_embeddings}")

    args.output_metrics.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(metrics).to_parquet(args.output_metrics, compression="zstd")
    print(f"Saved metrics -> {args.output_metrics}")

    if args.output_test_preds is not None:
        args.output_test_preds.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(test_preds_records).to_parquet(args.output_test_preds, compression="zstd")
        print(f"Saved test predictions -> {args.output_test_preds}")

    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
