"""Figure 5: Sankey diagram showing data flow from raw sources to ATLAS graph."""
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


def draw_sankey_segment(ax, x0, y0, x1, y1, width0, width1, color, alpha=0.6):
    """Draw a tapered polygon from (x0, y0±width0/2) to (x1, y1±width1/2)."""
    polygon = plt.Polygon([
        [x0, y0 - width0 / 2],
        [x0, y0 + width0 / 2],
        [x1, y1 + width1 / 2],
        [x1, y1 - width1 / 2],
    ], closed=True, facecolor=color, edgecolor="black", linewidth=0.5, alpha=alpha)
    ax.add_patch(polygon)


def main():
    outdir = Path("docs/figures")
    outdir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6)
    ax.axis("off")

    # Column positions
    col_x = [1, 4, 7, 10]
    labels = ["Raw Sources", "Processed", "Graph", "Outputs"]
    for x, lab in zip(col_x, labels):
        ax.text(x, 5.5, lab, ha="center", va="center", fontsize=12,
                fontweight="bold", color="dimgray")

    # Row positions and widths (relative)
    sources = {
        "UniProt": (1, 4.0, 0.8),
        "PDB": (1, 2.5, 0.5),
        "AlphaFold": (1, 1.5, 0.6),
        "Pfam": (1, 0.5, 0.4),
    }
    processed = {
        "targets_v2": (4, 4.0, 0.9),
        "structures": (4, 2.5, 0.7),
        "domains": (4, 1.2, 0.5),
    }
    graph_nodes = {
        "Protein": (7, 4.0, 0.9),
        "Structure": (7, 2.5, 0.6),
        "Domain": (7, 1.2, 0.5),
        "System": (7, 0.3, 0.3),
    }
    outputs = {
        "embeddings": (10, 4.0, 0.7),
        "atlas.gpickle": (10, 2.5, 0.6),
        "selection API": (10, 1.0, 0.5),
    }

    colors = plt.cm.Set3(np.linspace(0, 1, 10))

    # Draw nodes
    all_nodes = {**sources, **processed, **graph_nodes, **outputs}
    for i, (name, (x, y, w)) in enumerate(all_nodes.items()):
        rect = mpatches.FancyBboxPatch((x - 0.6, y - w / 2), 1.2, w,
                                       boxstyle="round,pad=0.02",
                                       facecolor=colors[i % len(colors)],
                                       edgecolor="black", linewidth=1)
        ax.add_patch(rect)
        ax.text(x, y, name, ha="center", va="center", fontsize=8, fontweight="bold")

    # Draw flows
    draw_sankey_segment(ax, 1.6, 4.0, 3.4, 4.0, 0.8, 0.9, colors[0])
    draw_sankey_segment(ax, 1.6, 2.5, 3.4, 2.5, 0.5, 0.7, colors[1])
    draw_sankey_segment(ax, 1.6, 1.5, 3.4, 2.5, 0.6, 0.7, colors[2])
    draw_sankey_segment(ax, 1.6, 0.5, 3.4, 1.2, 0.4, 0.5, colors[3])

    draw_sankey_segment(ax, 4.6, 4.0, 6.4, 4.0, 0.9, 0.9, colors[0])
    draw_sankey_segment(ax, 4.6, 2.5, 6.4, 2.5, 0.7, 0.6, colors[1])
    draw_sankey_segment(ax, 4.6, 1.2, 6.4, 1.2, 0.5, 0.5, colors[3])

    draw_sankey_segment(ax, 7.6, 4.0, 9.4, 4.0, 0.9, 0.7, colors[0])
    draw_sankey_segment(ax, 7.6, 2.5, 9.4, 2.5, 0.6, 0.6, colors[1])
    draw_sankey_segment(ax, 7.6, 0.3, 9.4, 1.0, 0.3, 0.5, colors[7])

    ax.set_title("ATLAS Data Flow: From Raw Sources to Knowledge Graph", fontsize=14,
                 fontweight="bold", pad=10)

    fig.savefig(outdir / "fig5_sankey.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "fig5_sankey.png", dpi=300, bbox_inches="tight")
    print(f"Saved {outdir}/fig5_sankey.pdf")


if __name__ == "__main__":
    main()
