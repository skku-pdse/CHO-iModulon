"""Approved Figure 5 plotting style; data are passed by the caller."""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
MAGENTA = "#CC79A7"
SKY = "#56B4E9"
DARK = "#242424"
GRAY = "#A7A7A7"
KEY_AXES = {
    87: {"short": "ISG-1", "color": BLUE},
    54: {"short": "ATF4", "color": ORANGE},
    39: {"short": "TDS", "color": GREEN},
    52: {"short": "ECM-1", "color": MAGENTA},
    65: {"short": "ECM-2", "color": SKY},
}


def set_publication_style() -> None:
    """Apply a compact style intended for two-column journal figures."""
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8.5,
            "axes.titlesize": 9.2,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.8,
            "ytick.labelsize": 7.8,
            "legend.fontsize": 7.5,
            "axes.linewidth": 0.7,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )


def clean_axis(ax: plt.Axes) -> None:
    """Use a restrained manuscript axis style."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out", pad=2)


def add_panel_letter(
    ax: plt.Axes, letter: str, x: float = -0.18, y: float = 1.08
) -> None:
    ax.text(
        x,
        y,
        letter,
        transform=ax.transAxes,
        fontsize=12.5,
        fontweight="bold",
        ha="left",
        va="top",
        clip_on=False,
    )


def plot_panel_a(ax: plt.Axes, best: pd.DataFrame, letter: bool = True) -> None:
    """Figure 5A: ranked best gene-weight match for every discovery iModulon."""
    ranked = best.sort_values("abs_pearson_r").reset_index(drop=True).copy()
    ranked["rank"] = np.arange(1, len(ranked) + 1)
    mnn = ranked["mutual_nearest_neighbor"].astype(bool)
    ax.scatter(
        ranked.loc[~mnn, "rank"],
        ranked.loc[~mnn, "abs_pearson_r"],
        s=15,
        facecolors="white",
        edgecolors=GRAY,
        linewidths=0.8,
        zorder=2,
    )
    ax.scatter(
        ranked.loc[mnn, "rank"],
        ranked.loc[mnn, "abs_pearson_r"],
        s=13,
        c=DARK,
        edgecolors="none",
        alpha=0.82,
        zorder=2,
    )
    offsets = {
        87: (-25, -14),
        54: (-12, -32),
        39: (-10, -17),
        52: (-25, 8),
        65: (-10, -17),
    }
    for discovery_id, info in KEY_AXES.items():
        row = ranked.loc[
            ranked["discovery_iModulon"].astype(int).eq(discovery_id)
        ].iloc[0]
        x = float(row["rank"])
        y = float(row["abs_pearson_r"])
        ax.scatter(
            [x],
            [y],
            s=28,
            c=info["color"],
            edgecolors="white",
            linewidths=0.55,
            zorder=4,
        )
        ax.annotate(
            info["short"],
            xy=(x, y),
            xytext=offsets[discovery_id],
            textcoords="offset points",
            color=info["color"],
            fontsize=8.7,
            fontweight="bold",
            arrowprops={
                "arrowstyle": "-",
                "color": info["color"],
                "lw": 1,
                "shrinkA": 1,
                "shrinkB": 1,
            },
            zorder=5,
        )
    ax.set_xlim(0, len(ranked) + 2)
    ax.set_ylim(0, 1.035)
    ax.set_xticks([1, 25, 50, 75, len(ranked)])
    ax.set_xlabel("Discovery iModulons ranked by best match")
    ax.set_ylabel("Best gene-weight |r|")
    clean_axis(ax)
    n_mnn = int(mnn.sum())
    ax.text(
        0.03,
        0.92,
        f"{n_mnn}/{len(ranked)} reciprocal\nbest matches",
        transform=ax.transAxes,
        ha="left",
        va="center",
        fontsize=7.6,
        color=DARK,
    )
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor=DARK,
            markeredgecolor=DARK,
            markersize=3.7,
            label="Reciprocal",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="white",
            markeredgecolor=GRAY,
            markersize=3.7,
            label="Non-reciprocal",
        ),
    ]
    ax.legend(
        handles=handles,
        loc="lower right",
        frameon=False,
        handletextpad=0.35,
        borderaxespad=0.25,
    )
    if letter:
        add_panel_letter(ax, "A")


plot_panel_b = plot_panel_a
