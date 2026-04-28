"""Figure 3: Primary benchmark bar chart with bootstrap CIs."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main():
    outdir = Path("docs/figures")
    outdir.mkdir(parents=True, exist_ok=True)

    # Primary benchmark: Protein→Domain link prediction, 20% hold-out test set.
    # GNN CIs: Mann-Whitney SE. Node2Vec CI: 1000× bootstrap on test set.
    # RBF/Quantum: node classification task (separate table in paper).
    models = ["GAT", "GraphSAGE", "Classical RBF*", "Quantum Kernel*", "Node2Vec"]
    auroc = [0.9705, 0.9664, 0.9331, 0.8731, 0.8411]
    ci_low = [0.9446, 0.9405, np.nan, np.nan, 0.8052]
    ci_high = [0.9964, 0.9923, np.nan, np.nan, 0.8342]
    colors = ["#2ecc71", "#3498db", "#95a5a6", "#9b59b6", "#e74c3c"]

    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(models))
    bars = ax.bar(x, auroc, color=colors, edgecolor="black", linewidth=1.2)

    # Error bars for CIs
    for i, (low, high) in enumerate(zip(ci_low, ci_high)):
        if not np.isnan(low):
            ax.plot([i, i], [low, high], "k-", linewidth=2)
            ax.plot([i], [low], "kv", markersize=6)
            ax.plot([i], [high], "k^", markersize=6)

    ax.set_ylabel("AUROC", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(models, rotation=15, ha="right", fontsize=10)
    ax.set_ylim(0.75, 1.0)
    ax.axhline(y=0.95, color="gray", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.text(4.5, 0.952, "0.95", fontsize=8, color="gray")

    for bar, val in zip(bars, auroc):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                f"{val:.4f}", ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax.set_title("Primary Benchmark: Protein→Domain Link Prediction\n* node classification task (Supplementary)", fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(outdir / "fig3_benchmark.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "fig3_benchmark.png", dpi=300, bbox_inches="tight")
    print(f"Saved {outdir}/fig3_benchmark.pdf")


if __name__ == "__main__":
    main()
