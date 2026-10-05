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
MID_GRAY = "#777777"
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


def plot_panel_b(ax: plt.Axes, activity: pd.DataFrame, letter: bool = True) -> None:
    """Figure 5B: gene-weight similarity versus activity conservation."""
    frame = activity.loc[activity["mutual_nearest_neighbor"].astype(bool)].copy()
    x = frame["M_abs_pearson_r"].to_numpy(dtype=float)
    y = frame["activity_pearson_r_sign_aligned"].to_numpy(dtype=float)
    relationship_r = float(np.corrcoef(x, y)[0, 1])
    ax.scatter(
        x, y, s=15, c=MID_GRAY, alpha=0.55, edgecolors="none", rasterized=True, zorder=2
    )
    offsets = {
        87: (-15, -20),
        54: (-50, 9),
        39: (-30, -30),
        52: (-30, 8),
        65: (-21, -30),
    }
    for discovery_id, info in KEY_AXES.items():
        row = frame.loc[frame["discovery_iModulon"].astype(int).eq(discovery_id)].iloc[
            0
        ]
        px = float(row["M_abs_pearson_r"])
        py = float(row["activity_pearson_r_sign_aligned"])
        ax.scatter(
            [px],
            [py],
            s=30,
            c=info["color"],
            edgecolors="white",
            linewidths=0.55,
            zorder=4,
        )
        if discovery_id not in offsets:
            continue
        ax.annotate(
            info["short"],
            xy=(px, py),
            xytext=offsets[discovery_id],
            textcoords="offset points",
            fontsize=8.7,
            fontweight="bold",
            color=info["color"],
            arrowprops={
                "arrowstyle": "-",
                "color": info["color"],
                "lw": 1,
                "shrinkA": 1,
                "shrinkB": 2,
            },
        )
    x_min = min(0.0, float(np.nanmin(x)) - 0.04)
    y_min = min(0.0, float(np.nanmin(y)) - 0.04)
    ax.set_xlim(x_min, 1.025)
    ax.set_ylim(y_min, 1.025)
    ax.set_xlabel("Gene-weight similarity, |r|")
    ax.set_ylabel("Sign-aligned activity Pearson r")
    clean_axis(ax)
    ax.text(
        0.97,
        0.08,
        f"Similarity–activity association\nacross reciprocal pairs: r = {relationship_r:.2f}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.6,
        color=DARK,
    )
    inset = ax.inset_axes([0.08, 0.6, 0.35, 0.41])
    overall = frame["activity_pearson_r_sign_aligned"].dropna().to_numpy()
    centered = frame["within_project_centered_pearson_r"].dropna().to_numpy()
    box = inset.boxplot(
        [overall, centered],
        positions=[1, 2],
        widths=0.52,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": DARK, "lw": 0.8},
        whiskerprops={"color": MID_GRAY, "lw": 0.6},
        capprops={"color": MID_GRAY, "lw": 0.6},
        boxprops={"edgecolor": MID_GRAY, "lw": 0.6},
    )
    box["boxes"][0].set_facecolor("#D8EAF5")
    box["boxes"][1].set_facecolor("#E5D9F2")
    rng = np.random.default_rng(41)
    for position, values, color in [(1, overall, BLUE), (2, centered, MAGENTA)]:
        jitter = rng.uniform(-0.12, 0.12, size=len(values))
        inset.scatter(
            position + jitter,
            values,
            s=3.2,
            c=color,
            alpha=0.28,
            edgecolors="none",
            rasterized=True,
        )
    inset.set_xticks([1, 2], ["Overall", "Project-\ncentered"])
    inset.set_ylim(y_min, 1.025)
    inset.set_title("Activity correlation", fontsize=7.0, pad=1.5)
    inset.tick_params(labelsize=6.3, length=1.7, pad=1.2)
    inset.spines["top"].set_visible(False)
    inset.spines["right"].set_visible(False)
    inset.spines["left"].set_linewidth(0.5)
    inset.spines["bottom"].set_linewidth(0.5)
    if letter:
        add_panel_letter(ax, "B", x=-0.16)


plot_panel_c = plot_panel_b
