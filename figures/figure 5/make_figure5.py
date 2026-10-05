"""Original A/B and lineage panels; newer member-colored C/D/E panels."""

from pathlib import Path
import ast
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import panel_a as pa
import panel_b as pb
import panel_f as pf
import panel_g as pg
from current_data import (
    load_matched_data,
    load_lineage_data,
    load_best_matches,
    load_activity,
)
from figure5_paths import RESULTS
from membership_panel import panelmatch

HERE = RESULTS
data = load_matched_data()
DM = data["discovery_m"]
EM = data["expanded_m"].loc[DM.index]
dg = data["discovery_members"]
eg = data["expanded_members"]
t = pd.read_csv(RESULTS / "activity_and_membership_comparison.csv")


def draw_membership(ax, d, label):
    panelmatch(ax, d, label, DM, EM, dg, eg, t)


def legend(fig):
    fig.legend(
        handles=[
            Line2D([], [], marker="o", ls="", color=c, label=n)
            for c, n in [
                ("#0072B2", "Shared member"),
                ("#E69F00", "Discovery only"),
                ("#009E73", "Expanded only"),
            ]
        ]
        + [
            Line2D(
                [],
                [],
                marker="D",
                ls="",
                markerfacecolor="none",
                color="#6F4AA8",
                label="Four histone genes",
            )
        ],
        loc="lower center",
        bbox_to_anchor=(0.53, 0.005),
        ncol=4,
        frameon=False,
        fontsize=7,
    )


def save(fig, stem):
    fig.canvas.draw()
    for ax in fig.axes:
        if ax.get_box_aspect() == 1:
            bounds = ax.get_window_extent()
            assert abs(bounds.width / bounds.height - 1) < 1e-06, "Panel must be square"
    for ext in ["png", "pdf"]:
        fig.savefig(
            HERE / f"{stem}.{ext}", dpi=400, bbox_inches="tight", pad_inches=0.05
        )
    plt.close(fig)


pa.set_publication_style()
for d, label in [(87, "C"), (54, "D"), (39, "E")]:
    fig, ax = plt.subplots(figsize=(3.8, 3.8), layout="constrained")
    draw_membership(ax, d, label)
    ax.set_box_aspect(1)
    fig.legend(
        handles=[
            Line2D([], [], marker="o", ls="", color=c, label=n)
            for c, n in [
                ("#0072B2", "Shared"),
                ("#E69F00", "Discovery only"),
                ("#009E73", "Expanded only"),
            ]
        ],
        loc="outside lower center",
        ncol=3,
        frameon=False,
        fontsize=7,
    )
    save(fig, f"figure6{label}_membership_hybrid")
fig = plt.figure(figsize=(8.1, 9.7))
outer = fig.add_gridspec(
    3,
    1,
    height_ratios=[1.02, 1.18, 0.9],
    left=0.09,
    right=0.985,
    bottom=0.08,
    top=0.97,
    hspace=0.52,
)
top = outer[0, 0].subgridspec(1, 2, width_ratios=[0.92, 1.08], wspace=0.34)
pa.plot_panel_a(fig.add_subplot(top[0, 0]), load_best_matches())
pb.plot_panel_b(fig.add_subplot(top[0, 1]), load_activity())
middle = outer[1, 0].subgridspec(1, 3, wspace=0.62)
for j, (d, label) in enumerate([(87, "C"), (54, "D"), (39, "E")]):
    ax = fig.add_subplot(middle[0, j])
    draw_membership(ax, d, label)
    ax.set_box_aspect(1)
bottom = outer[2, 0].subgridspec(1, 2, width_ratios=[2.45, 1.0], wspace=0.4)
axf = fig.add_subplot(bottom[0, 0])
axg = fig.add_subplot(bottom[0, 1])
ld = load_lineage_data()
pf.plot_im86_activity(axf, ld)
pf.add_panel_letter(axf, "F", x=-0.16, y=1.1)
pg.plot_panel_e(axg, ld, letter=False)
pg.add_panel_letter(axg, "G")
legend(fig)
save(fig, "figure5_compendium_comparison")
print("Saved hybrid Figure 5 and member-colored C/D/E panels.")
