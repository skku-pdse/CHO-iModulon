"""Approved Figure 5 plotting style; data are passed by the caller."""

from __future__ import annotations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

WITHIN_PROJECT = "PRJNA1130622"
ORANGE = "#E69F00"
ECM4_COLOR = "#6F4AA8"
DARK = "#242424"
GRAY = "#A7A7A7"
PALE_GRAY = "#EEEEEE"


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
            "legend.fontsize": 8,
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


def plot_panel_e(ax: plt.Axes, data: dict, letter: bool = True) -> None:
    """Figure 5E: within-project K1 versus DXB11 comparison for PRJNA1130622."""
    rng = np.random.default_rng(91)
    ax.axhline(0, color=PALE_GRAY, lw=0.55, zorder=0)
    within = (
        data["im86_samples"]
        .loc[
            data["im86_samples"]["Project"].astype(str).eq(WITHIN_PROJECT)
            & data["im86_samples"]["Cell line"].isin(["CHO-K1", "CHO-DXB11"])
        ]
        .copy()
    )
    group_sizes = {}
    for position, cell_line in enumerate(["CHO-K1", "CHO-DXB11"], start=1):
        values = (
            within.loc[within["Cell line"].eq(cell_line), "iM86_activity"]
            .dropna()
            .to_numpy()
        )
        group_sizes[cell_line] = len(values)
        jitter = rng.uniform(-0.11, 0.11, size=len(values))
        ax.scatter(
            position + jitter,
            values,
            s=14,
            c=GRAY if cell_line == "CHO-K1" else ECM4_COLOR,
            alpha=0.78,
            edgecolors="white",
            linewidths=0.3,
            rasterized=True,
        )
        ax.plot(
            [position - 0.2, position + 0.2],
            [np.median(values), np.median(values)],
            color=DARK,
            lw=1.2,
        )
    ax.set_xlim(0.5, 2.5)
    ax.set_xticks(
        [1, 2],
        [
            f"CHO-K1\n(n={group_sizes['CHO-K1']})",
            f"CHO-DXB11\n(n={group_sizes['CHO-DXB11']})",
        ],
    )
    ax.set_ylim(-15, 13)
    ax.set_ylabel("Expanded iM86 activity")
    ax.set_title(WITHIN_PROJECT, loc="center", pad=3)
    ax.text(
        0.5,
        0.5,
        "Exact permutation\n" + format_p_value(data["permutation_p"]),
        transform=ax.transAxes,
        ha="center",
        va="center",
        color=DARK,
    )
    clean_axis(ax)
    if letter:
        add_panel_letter(ax, "E")


plot_panel_f = plot_panel_e
