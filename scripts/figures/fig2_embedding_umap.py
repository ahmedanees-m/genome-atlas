"""Figure 2: UMAP visualization of protein embeddings colored by system type."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    outdir = Path("docs/figures")
    outdir.mkdir(parents=True, exist_ok=True)

    # Try to load real embeddings; fall back to synthetic demo
    emb_path = Path.home() / "pen-stack/data/embeddings/graphsage.parquet"
    if emb_path.exists():
        df = pd.read_parquet(emb_path)
        print(f"Loaded {len(df)} embeddings from {emb_path}")
    else:
        print("Real embeddings not found — generating synthetic demo")
        np.random.seed(42)
        n = 500
        df = pd.DataFrame({
            "node_id": [f"Protein_demo_{i}" for i in range(n)],
            "embedding": [np.random.randn(128).astype(np.float32) for _ in range(n)],
            "system_type": np.random.choice(
                ["CRISPR-Cas", "Tyrosine_recombinase", "Prime_editor",
                 "CRISPR-Transposon", "DDE_transposase", "Bridge_Recombinase"],
                size=n
            ),
        })

    # UMAP reduce
    try:
        import umap
        reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
        embeddings = np.vstack(df["embedding"].values)
        xy = reducer.fit_transform(embeddings)
    except ImportError:
        print("umap-learn not installed; using PCA fallback")
        from sklearn.decomposition import PCA
        embeddings = np.vstack(df["embedding"].values)
        xy = PCA(n_components=2).fit_transform(embeddings)

    fig, ax = plt.subplots(figsize=(8, 8))
    cmap = {
        "CRISPR-Cas": "#e74c3c",
        "Tyrosine_recombinase": "#3498db",
        "Prime_editor": "#2ecc71",
        "CRISPR-Transposon": "#9b59b6",
        "DDE_transposase": "#f39c12",
        "Bridge_Recombinase": "#1abc9c",
    }
    for stype, color in cmap.items():
        mask = df["system_type"] == stype
        if mask.sum() == 0:
            continue
        ax.scatter(xy[mask, 0], xy[mask, 1], c=color, label=stype, s=8, alpha=0.6)

    ax.set_xlabel("UMAP 1", fontsize=12)
    ax.set_ylabel("UMAP 2", fontsize=12)
    ax.set_title("ESM-2 + GraphSAGE Embedding Space (colored by system type)",
                 fontsize=14, fontweight="bold")
    ax.legend(loc="best", fontsize=8, frameon=True)
    fig.savefig(outdir / "fig2_umap.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "fig2_umap.png", dpi=300, bbox_inches="tight")
    print(f"Saved {outdir}/fig2_umap.pdf")


if __name__ == "__main__":
    main()
