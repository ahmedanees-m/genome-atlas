"""Figure 4: Radar chart comparing selection axes for three example systems."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main():
    outdir = Path("docs/figures")
    outdir.mkdir(parents=True, exist_ok=True)

    categories = ["DSB Avoidance", "AAV Fit", "Cargo Size", "Cell Type"]
    N = len(categories)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]

    systems = {
        "SpCas9": [0.2, 0.1, 0.9, 0.95],
        "PE2": [0.85, 0.1, 1.0, 0.9],
        "CAST-I-F_evoCAST": [1.0, 0.1, 1.0, 0.85],
        "Cas12f": [0.2, 1.0, 0.9, 0.95],
    }
    colors = ["#e74c3c", "#2ecc71", "#9b59b6", "#3498db"]

    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))

    for (name, values), color in zip(systems.items(), colors):
        values += values[:1]
        ax.plot(angles, values, "o-", linewidth=2, label=name, color=color)
        ax.fill(angles, values, alpha=0.15, color=color)

    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories, fontsize=10)
    ax.set_ylim(0, 1)
    ax.set_title("Selection Axes Comparison (Example Systems)", fontsize=13,
                 fontweight="bold", pad=20)
    ax.legend(loc="upper right", bbox_to_anchor=(1.3, 1.1), fontsize=9)

    fig.savefig(outdir / "fig4_radar.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(outdir / "fig4_radar.png", dpi=300, bbox_inches="tight")
    print(f"Saved {outdir}/fig4_radar.pdf")


if __name__ == "__main__":
    main()
