"""Figure 1: ATLAS knowledge graph architecture schema."""
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


def draw_node(ax, x, y, label, color, size=0.08):
    circle = plt.Circle((x, y), size, color=color, ec="black", linewidth=1.5, zorder=3)
    ax.add_patch(circle)
    ax.text(x, y, label, ha="center", va="center", fontsize=8, fontweight="bold", zorder=4)


def draw_edge(ax, x1, y1, x2, y2, label, color="gray"):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="->", color=color, lw=1.5))
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    ax.text(mx, my + 0.05, label, ha="center", va="bottom", fontsize=7,
            style="italic", color=color)


def main():
    outdir = Path("docs/figures")
    outdir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")

    # Node positions
    nodes = {
        "System": (2, 5),
        "Protein": (5, 5),
        "Domain": (8, 5),
        "Structure": (5, 2.5),
        "Organism": (2, 2.5),
    }
    colors = {
        "System": "#e74c3c",
        "Protein": "#3498db",
        "Domain": "#2ecc71",
        "Structure": "#9b59b6",
        "Organism": "#f39c12",
    }

    for name, (x, y) in nodes.items():
        draw_node(ax, x, y, name, colors[name])

    # Edges
    draw_edge(ax, 2, 4.8, 4.5, 4.8, "HAS_PROTEIN")
    draw_edge(ax, 5.5, 4.8, 7.5, 4.8, "HAS_DOMAIN")
    draw_edge(ax, 5, 4.5, 5, 3.0, "HAS_STRUCTURE")
    draw_edge(ax, 2.3, 4.5, 2, 3.0, "DERIVED_FROM")
    draw_edge(ax, 5, 4.5, 2.3, 3.0, "MECHANISM_OF")

    # Counts
    counts = {
        "System": "16",
        "Protein": "9,500",
        "Domain": "18 Pfams",
        "Structure": "2,239",
        "Organism": "347",
    }
    for name, (x, y) in nodes.items():
        ax.text(x, y - 0.25, counts[name], ha="center", va="top", fontsize=8, color="dimgray")

    ax.set_title("GENOME-ATLAS Heterogeneous Knowledge Graph Schema (v0.5.1)",
                 fontsize=14, fontweight="bold", pad=20)

    # Legend
    legend_patches = [mpatches.Patch(color=c, label=n) for n, c in colors.items()]
    ax.legend(handles=legend_patches, loc="lower right", frameon=True, fontsize=8)

    fig.savefig(outdir / "fig1_architecture.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "fig1_architecture.png", dpi=300, bbox_inches="tight")
    print(f"Saved {outdir}/fig1_architecture.pdf")


if __name__ == "__main__":
    main()
