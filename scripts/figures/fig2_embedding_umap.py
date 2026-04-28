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
        # Subset to Protein nodes only for cleaner visualization (9,500 points)
        if "node_type" in df.columns:
            df = df[df["node_type"] == "Protein"].copy()
        # Color by mechanism_bucket via targets_v2 if available
        targets_path = Path.home() / "pen-stack/data/processed/targets_v2.parquet"
        label_col = "node_type"
        if targets_path.exists() and "node_id" in df.columns:
            targets = pd.read_parquet(targets_path, columns=["accession", "primary_mechanism_bucket"])
            targets["node_id"] = "Protein_" + targets["accession"].astype(str)
            df = df.merge(targets[["node_id", "primary_mechanism_bucket"]],
                          on="node_id", how="left")
            df["primary_mechanism_bucket"] = df["primary_mechanism_bucket"].fillna("Unknown")
            label_col = "primary_mechanism_bucket"
    else:
        print("Real embeddings not found — generating synthetic demo")
        np.random.seed(42)
        n = 500
        df = pd.DataFrame({
            "node_id": [f"Protein_demo_{i}" for i in range(n)],
            "embedding": [np.random.randn(128).astype(np.float32) for _ in range(n)],
            "primary_mechanism_bucket": np.random.choice(
                ["DSB_NUCLEASE", "DSB_FREE_TRANSEST_RECOMBINASE", "TRANSPOSASE"],
                size=n
            ),
        })
        label_col = "primary_mechanism_bucket"

    # UMAP reduce (fall back to PCA if umap-learn not installed)
    try:
        import umap
        reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
        embeddings = np.vstack(df["embedding"].values)
        xy = reducer.fit_transform(embeddings)
        dim_label = "UMAP"
    except ImportError:
        print("umap-learn not installed; using PCA fallback")
        from sklearn.decomposition import PCA
        embeddings = np.vstack(df["embedding"].values)
        xy = PCA(n_components=2).fit_transform(embeddings)
        dim_label = "PC"

    fig, ax = plt.subplots(figsize=(8, 8))
    cmap = {
        "DSB_NUCLEASE": "#e74c3c",
        "DSB_FREE_TRANSEST_RECOMBINASE": "#3498db",
        "TRANSPOSASE": "#f39c12",
        "Unknown": "#cccccc",
        # Legacy system_type values (demo)
        "CRISPR-Cas": "#e74c3c",
        "Tyrosine_recombinase": "#3498db",
        "Prime_editor": "#2ecc71",
        "CRISPR-Transposon": "#9b59b6",
        "DDE_transposase": "#f39c12",
        "Bridge_Recombinase": "#1abc9c",
    }
    for label in df[label_col].unique():
        color = cmap.get(label, "#7f8c8d")
        mask = df[label_col] == label
        ax.scatter(xy[mask, 0], xy[mask, 1], c=color, label=label, s=4, alpha=0.5)

    ax.set_xlabel(f"{dim_label} 1", fontsize=12)
    ax.set_ylabel(f"{dim_label} 2", fontsize=12)
    ax.set_title("GraphSAGE Embedding Space (proteins colored by mechanism bucket)",
                 fontsize=13, fontweight="bold")
    ax.legend(loc="best", fontsize=9, frameon=True, markerscale=3)
    fig.savefig(outdir / "fig2_umap.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "fig2_umap.png", dpi=300, bbox_inches="tight")
    print(f"Saved {outdir}/fig2_umap.pdf")


if __name__ == "__main__":
    main()
