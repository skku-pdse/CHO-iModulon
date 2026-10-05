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
from figure2_paths import DATA, RESULTS

TARGET_IM = 52
TARGET_IM_NAME = "Lineage / ECM-adhesion"
POS_COLOR = "#3b82bd"
NEG_COLOR = "#c8553d"
FIG_W, FIG_H = 4.6, 4.2
DPI = 300
BAR_HEIGHT = 0.72
LABEL_FONT = 8.0
VALUE_FONT = 6.8
TICK_FONT = 7.5


def plot_panel(genes, n_total, out_stem):
    g = genes.sort_values("Loading", ascending=True).reset_index(drop=True)
    n = len(g)
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    y = np.arange(n)
    colors = [POS_COLOR if v >= 0 else NEG_COLOR for v in g["Loading"].values]
    ax.barh(
        y,
        g["Loading"].values,
        height=BAR_HEIGHT,
        color=colors,
        edgecolor="white",
        linewidth=0.5,
    )
    ax.axvline(0, color="0.3", linewidth=0.8, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels(g["DisplayLabel"].tolist(), fontsize=LABEL_FONT)
    for tick_label, label_disp in zip(ax.get_yticklabels(), g["DisplayLabel"].values):
        is_loc_fallback = str(label_disp).startswith("LOC")
        if not is_loc_fallback:
            tick_label.set_fontstyle("italic")
    xmin, xmax = (g["Loading"].min(), g["Loading"].max())
    pad = (xmax - xmin) * 0.015
    for yi, v in zip(y, g["Loading"].values):
        ha = "left" if v >= 0 else "right"
        x_text = v + pad if v >= 0 else v - pad
        ax.text(
            x_text,
            yi,
            f"{v:+.3f}",
            va="center",
            ha=ha,
            fontsize=VALUE_FONT,
            color="0.25",
        )
    ax.set_xlabel("Gene loading (M coefficient)", fontsize=9)
    title_n = f"top {n} of {n_total}" if n < n_total else f"all {n}"
    ax.set_title(
        f"iM{TARGET_IM}: {TARGET_IM_NAME}\n{title_n} gene members by |loading|",
        fontsize=10,
        pad=8,
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="x", labelsize=TICK_FONT)
    ax.grid(axis="x", linestyle=":", alpha=0.4)
    span = xmax - xmin
    ax.set_xlim(xmin - span * 0.18, xmax + span * 0.18)
    ax.set_ylim(-0.7, n - 0.3)
    handles = [
        Patch(facecolor=POS_COLOR, edgecolor="white", label="positive loading"),
        Patch(facecolor=NEG_COLOR, edgecolor="white", label="negative loading"),
    ]
    ax.legend(
        handles=handles,
        loc="lower right",
        fontsize=7.0,
        frameon=True,
        framealpha=0.95,
        borderpad=0.4,
        handletextpad=0.5,
        handlelength=1.2,
    )
    fig.tight_layout()
    fig.savefig(out_stem.with_suffix(".png"), dpi=DPI, bbox_inches="tight")
    fig.savefig(out_stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    genes = pd.read_csv(DATA / "iM52_figure_genes.csv")
    output = RESULTS / "gene_loadings"
    output.mkdir(exist_ok=True)
    plot_panel(genes, 38, output / "iM52_topgenes")
