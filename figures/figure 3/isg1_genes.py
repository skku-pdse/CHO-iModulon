"""Plot frozen manuscript gene labels; no annotation or symbol matching runs here.
Starred labels retain their original curated display meaning and are not new
evidence of sequence identity. Gene identifiers remain in the input table.
"""

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from figure3_paths import DATA, RESULTS

TARGET_IM = 87
TARGET_IM_NAME = "Antiviral / type-I IFN response"
BAR_COLOR = "#3b82bd"
CANONICAL_MARKERS = {}
FIG_W, FIG_H = 4.6, 3.9
DPI = 300
BAR_HEIGHT = 0.72
LABEL_FONT = 8.0
VALUE_FONT = 6.8
TICK_FONT = 7.5


def plot_panel(genes, out_stem):
    g = genes.sort_values("Loading", ascending=True).reset_index(drop=True)
    n = len(g)
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    y = np.arange(n)
    ax.barh(
        y,
        g["Loading"].values,
        height=BAR_HEIGHT,
        color=BAR_COLOR,
        edgecolor="white",
        linewidth=0.5,
    )
    ax.axvline(0, color="0.3", linewidth=0.8, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(g["DisplayLabel"].tolist(), fontsize=LABEL_FONT)
    for tick_label, gene_orig, label_disp in zip(
        ax.get_yticklabels(), g["Gene"].values, g["DisplayLabel"].values
    ):
        is_loc_fallback = str(label_disp).startswith("LOC")
        if not is_loc_fallback:
            tick_label.set_fontstyle("italic")
        if (
            str(gene_orig) in CANONICAL_MARKERS
            or str(label_disp).rstrip("*") in CANONICAL_MARKERS
        ):
            tick_label.set_fontweight("bold")
            tick_label.set_color("black")
    xmin, xmax = (0.0, g["Loading"].max())
    pad = (xmax - xmin) * 0.012
    for yi, v in zip(y, g["Loading"].values):
        ax.text(
            v + pad,
            yi,
            f"{v:+.3f}",
            va="center",
            ha="left",
            fontsize=VALUE_FONT,
            color="0.25",
        )
    ax.set_xlabel("Gene loading (M coefficient)", fontsize=9)
    ax.set_title(
        f"Interferon-Stimulated Genes (ISG)-1 (iM{TARGET_IM})",
        fontsize=10,
        pad=8,
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=TICK_FONT)
    ax.grid(axis="x", linestyle=":", alpha=0.4)
    span = xmax - xmin
    ax.set_xlim(xmin - span * 0.02, xmax + span * 0.18)
    ax.set_ylim(-0.7, n - 0.3)
    fig.tight_layout()
    fig.savefig(out_stem.with_suffix(".png"), dpi=DPI, bbox_inches="tight")
    fig.savefig(out_stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    genes = pd.read_csv(DATA / "iM87_figure_genes.csv")
    output = RESULTS / "gene_loadings"
    output.mkdir(exist_ok=True)
    plot_panel(genes, output / "figure3F_ISG1_gene_loadings")
