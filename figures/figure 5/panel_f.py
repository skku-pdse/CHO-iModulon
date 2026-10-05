"""Approved Figure 5 plotting style; data are passed by the caller."""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
DARK = "#242424"
MID_GRAY = "#777777"
GRAY = "#A7A7A7"
PALE_GRAY = "#EEEEEE"
ECM4_COLOR = "#6F4AA8"
ECM_GENES = {"Efemp1", "Col3a1", "Flna", "Thbs1", "Pdpn", "Procr"}


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


def format_p_value(value: float) -> str:
    if value < 0.001:
        exponent = int(np.floor(np.log10(value)))
        coefficient = value / 10**exponent
        return f"$P={coefficient:.2f}\\times10^{{{exponent}}}$"
    return f"$P={value:.3f}$"


def plot_im86_activity(ax: plt.Axes, data: dict) -> None:
    frame = data["im86_samples"].copy()
    categories = [
        ("CHO-K1", "K1"),
        ("CHO-S", "S"),
        ("CHO-DG44", "DG44"),
        ("CHO-ZeLa", "ZeLa"),
        ("NIST-CHO", "NIST-\nCHO"),
        ("CHO-DXB11", "DXB11"),
    ]
    rng = np.random.default_rng(73)
    positions = np.arange(1, len(categories) + 1)
    activity_groups: list[np.ndarray] = []
    for position, (cell_line, _) in zip(positions, categories):
        values = frame.loc[
            frame["Cell line"].astype(str).eq(cell_line), "iM86_activity"
        ].dropna()
        activity_groups.append(values.to_numpy(dtype=float))
        if cell_line != "CHO-DXB11":
            jitter = rng.uniform(-0.18, 0.18, size=len(values))
            ax.scatter(
                position + jitter,
                values,
                s=4.0,
                c=GRAY,
                alpha=0.28,
                edgecolors="none",
                rasterized=True,
                zorder=1,
            )
        else:
            dxb = frame.loc[frame["Cell line"].astype(str).eq("CHO-DXB11")].copy()
            project_styles = [
                ("PRJEB38542", ECM4_COLOR, "o"),
                ("PRJNA1130622", ECM4_COLOR, "^"),
            ]
            for project, color, marker in project_styles:
                values_project = dxb.loc[
                    dxb["Project"].astype(str).eq(project), "iM86_activity"
                ].dropna()
                jitter = rng.uniform(-0.17, 0.17, size=len(values_project))
                ax.scatter(
                    position + jitter,
                    values_project,
                    s=12,
                    c=color,
                    marker=marker,
                    alpha=0.72,
                    edgecolors="white",
                    linewidths=0.25,
                    rasterized=True,
                    zorder=3,
                )
    box = ax.boxplot(
        activity_groups,
        positions=positions,
        widths=0.52,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": DARK, "lw": 0.8},
        whiskerprops={"color": MID_GRAY, "lw": 0.6},
        capprops={"color": MID_GRAY, "lw": 0.6},
        boxprops={"edgecolor": MID_GRAY, "lw": 0.6},
    )
    for index, artist in enumerate(box["boxes"]):
        artist.set_facecolor("#E8E8E8" if index < 5 else "#EDE5F5")
        artist.set_alpha(0.75)
    full_labels = [
        f"{cell_line}\n(n={len(values)})"
        for (cell_line, _), values in zip(categories, activity_groups)
    ]
    ax.set_xticks(positions, full_labels)
    plt.setp(ax.get_xticklabels(), rotation=0, ha="center", fontsize=7.2)
    ax.set_ylabel("Activity coefficient")
    ax.set_xlabel("")
    ax.set_title("Extracellular Region-8 (expanded iM86)", loc="center", pad=4)
    ax.set_ylim(-15, 13)
    clean_axis(ax)
    ax.legend(
        handles=[
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markerfacecolor=ECM4_COLOR,
                markeredgecolor="white",
                markersize=4.2,
                label="DXB11: PRJEB38542",
            ),
            Line2D(
                [0],
                [0],
                marker="^",
                linestyle="none",
                markerfacecolor=ECM4_COLOR,
                markeredgecolor="white",
                markersize=4.5,
                label="DXB11: PRJNA1130622",
            ),
        ],
        loc="lower left",
        frameon=False,
        ncol=1,
        handletextpad=0.3,
        borderaxespad=0.2,
        fontsize=7.5,
    )
