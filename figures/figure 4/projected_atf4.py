"""ATF4 projected activities using the original two-panel figure style.
Activities and metadata are joined by SampleID from the current validation input.
Original displayed samples agree with the historical CSV within its rounding.
"""

from figure4_paths import DATA, RESULTS
import matplotlib

matplotlib.use("Agg")
OUT_DIR = RESULTS / "external_projection"
OUT_DIR.mkdir(exist_ok=True)
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

TARGET_IM = 54
TARGET_NAME = "AAR"
COL_CTRL = "#2C9F8E"
COL_PERT = "#E67E22"
COL_BOX = "#2b2b2b"
FIG_W, FIG_H = (13.0, 4.6)
DPI = 300
DOT_SIZE = 22
MEAN_DOT_SIZE = 100
MEAN_EDGE_W = 1.4
DOT_ALPHA = 0.85
DOT_EDGEW = 0.5
LABEL_FONT = 14
TICK_FONT = 12
TITLE_FONT = 12
SUPTITLE_FONT = 11.5
ANNOT_FONT = 10


def draw_prjeb37009(ax, activity, m37):
    """
    4 columns: paired 8mM/0mM within parental and sub-clone groups.
        col 0  parental  8 mM Q (control)         filled circle
        col 1  parental  0 mM Q (perturbation)    filled circle
        --- visual gap ---
        col 2  sub-clone 8 mM Q (control)         open circle
        col 3  sub-clone 0 mM Q (perturbation)    open circle
    Hy_* samples (different cell-line background) are excluded.
    No box, no whiskers, no connecting line — mean shown as larger circle.
    """
    cond_col = "Condition "

    def pick(cond_list):
        out = []
        for cond in cond_list:
            sub = m37[m37[cond_col] == cond]
            if len(sub) == 0:
                raise RuntimeError(f"PRJEB37009: no sample for '{cond}'")
            for _, r in sub.iterrows():
                out.append(float(activity[r["col"]]))
        return out

    ctrl_parental = pick(["K1_8mMQ_R1", "K1_8mMQ_R2"])
    pert_parental = pick(["K1_0mMQ_R1", "K1_0mMQ_R2"])
    ctrl_sub2 = pick(["K1_8mMQ_sub2_R1", "K1_8mMQ_sub2_R2"])
    pert_sub2 = pick(["K1_0mMQ_sub2_R1", "K1_0mMQ_sub2_R2"])
    hy_n = m37[cond_col].astype(str).str.contains("Hy").sum()
    GROUP_GAP = 0.7
    x_par_8 = 0.0
    x_par_0 = 1.5
    x_sub_8 = 1.5 + 1.0 + GROUP_GAP
    x_sub_0 = x_sub_8 + 1.5
    columns = [
        (x_par_8, ctrl_parental, COL_CTRL, False),
        (x_par_0, pert_parental, COL_PERT, False),
        (x_sub_8, ctrl_sub2, COL_CTRL, True),
        (x_sub_0, pert_sub2, COL_PERT, True),
    ]
    rng = np.random.default_rng(54)
    for x, vals, color, is_sub in columns:
        vals = np.asarray(vals)
        if len(vals) == 0:
            continue
        for v in vals:
            xj = x + rng.uniform(-0.1, 0.1)
            ax.scatter(
                xj,
                v,
                s=DOT_SIZE,
                color=color,
                alpha=DOT_ALPHA,
                edgecolor="white",
                linewidth=DOT_EDGEW,
                marker="o",
                zorder=4,
            )
        m = float(vals.mean())
        if color == COL_PERT:
            ax.hlines(
                y=m, xmin=x - 0.75, xmax=x - 0.22, color=color, linewidth=2.2, zorder=5
            )
        else:
            ax.hlines(
                y=m, xmin=x + 0.22, xmax=x + 0.75, color=color, linewidth=2.2, zorder=5
            )
    ax.axhline(0, color="#666", linewidth=0.7, linestyle="--", alpha=0.6, zorder=1)
    xs = [x_par_8, x_par_0, x_sub_8, x_sub_0]
    group_mids = [(x_par_8 + x_par_0) / 2, (x_sub_8 + x_sub_0) / 2]
    ax.set_xticks(group_mids)
    ax.set_xticklabels(["parental", "sub-clone"], fontsize=TICK_FONT + 0.5)
    left_columns_for_labels = [
        (x_par_8, False, len(ctrl_parental)),
        (x_par_0, True, len(pert_parental)),
        (x_sub_8, False, len(ctrl_sub2)),
        (x_sub_0, True, len(pert_sub2)),
    ]
    yb_offset = -0.02
    for x, is_pert, n in left_columns_for_labels:
        q_disp = "8 mM Q" if not is_pert else "0 mM Q"
        col_color = COL_CTRL if not is_pert else COL_PERT
        ax.annotate(
            f"{q_disp}\nn={n}",
            xy=(x, yb_offset),
            xycoords=("data", "axes fraction"),
            xytext=(0, -16),
            textcoords="offset points",
            ha="center",
            va="top",
            fontsize=TICK_FONT - 0.5,
            color=col_color,
            annotation_clip=False,
        )
    ax.set_xlim(min(xs) - 0.55, max(xs) + 0.55)
    ax.set_ylabel(f"iM54 (ATF4) activity (A')", fontsize=LABEL_FONT)
    ax.set_title("PRJEB37009", fontsize=TITLE_FONT, pad=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="y", labelsize=TICK_FONT)


def draw_prjeb55916(ax, activity, m55):
    """
    6 paired columns: (6h_8mMQ, 6h_0mMQ) (24h_8mMQ, 24h_0mMQ) (48h_8mMQ, 48h_0mMQ)
    Boxplot baseline + ERR-level dots.
    """
    cond_col = "Condition "
    timepoints = [("6 h", "6h"), ("24 h", "24h"), ("48 h", "48h")]
    columns = []
    for tp_label, tp_key in timepoints:
        for is_pert, q_key in [(False, "8mMQ"), (True, "0mMQ")]:
            mask = (
                m55[cond_col]
                .astype(str)
                .str.contains(f"_{q_key}_{tp_key}_", regex=False)
            )
            sub = m55[mask]
            vals = sub["col"].apply(lambda c: float(activity[c])).tolist()
            columns.append(
                dict(tp_label=tp_label, q_key=q_key, is_pert=is_pert, vals=vals)
            )
    pair_gap = 1.0
    group_gap = 1.0
    xs = []
    for ti in range(3):
        x0 = ti * group_gap + ti * pair_gap
        xs.append(x0 + 0.0)
        xs.append(x0 + 1.0)
    xs = []
    cursor = 0.0
    for ti in range(3):
        xs.append(cursor)
        xs.append(cursor + 1.0)
        cursor += 1.0 + group_gap
    rng = np.random.default_rng(11)
    for x, c in zip(xs, columns):
        vals = np.asarray(c["vals"])
        if len(vals) == 0:
            continue
        color = COL_PERT if c["is_pert"] else COL_CTRL
        for v in vals:
            xj = x + rng.uniform(-0.16, 0.16)
            ax.scatter(
                xj,
                v,
                s=DOT_SIZE,
                color=color,
                alpha=DOT_ALPHA,
                edgecolor="white",
                linewidth=DOT_EDGEW,
                marker="o",
                zorder=4,
            )
        m = float(vals.mean())
        if c["is_pert"]:
            ax.hlines(
                y=m, xmin=x - 0.5, xmax=x - 0.22, color=color, linewidth=2.2, zorder=5
            )
        else:
            ax.hlines(
                y=m, xmin=x + 0.22, xmax=x + 0.5, color=color, linewidth=2.2, zorder=5
            )
    ax.axhline(0, color="#666", linewidth=0.7, linestyle="--", alpha=0.6, zorder=1)
    pair_mids = [(xs[2 * i] + xs[2 * i + 1]) / 2 for i in range(3)]
    ax.set_xticks(pair_mids)
    ax.set_xticklabels([tp[0] for tp in timepoints], fontsize=TICK_FONT + 0.5)
    yb_offset = -0.02
    for x, c in zip(xs, columns):
        q_disp = "8 mM Q" if not c["is_pert"] else "0 mM Q"
        col_color = COL_CTRL if not c["is_pert"] else COL_PERT
        n = len(c["vals"])
        ax.annotate(
            f"{q_disp}\nn={n}",
            xy=(x, yb_offset),
            xycoords=("data", "axes fraction"),
            xytext=(0, -16),
            textcoords="offset points",
            ha="center",
            va="top",
            fontsize=TICK_FONT - 0.5,
            color=col_color,
            annotation_clip=False,
        )
    ax.set_xlim(xs[0] - 0.7, xs[-1] + 0.7)
    ax.set_title("PRJEB55916", fontsize=TITLE_FONT, pad=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="y", labelsize=TICK_FONT)


def plot(A, df_meta):
    m37 = project_meta(A, df_meta, "PRJEB37009")
    m55 = project_meta(A, df_meta, "PRJEB55916")
    activity = A.loc[TARGET_IM]
    fig = plt.figure(figsize=(FIG_W, FIG_H))
    gs = GridSpec(1, 2, width_ratios=[0.85, 1.55], wspace=0.1, figure=fig)
    ax_left = fig.add_subplot(gs[0, 0])
    ax_right = fig.add_subplot(gs[0, 1], sharey=ax_left)
    focus_vals = []
    for c in m55["col"]:
        focus_vals.append(float(activity[c]))
    cond_col = "Condition "
    k1_mask = ~m37[cond_col].astype(str).str.contains("Hy")
    for c in m37[k1_mask]["col"]:
        focus_vals.append(float(activity[c]))
    focus_vals = np.array(focus_vals)
    ymin, ymax = (float(focus_vals.min()), float(focus_vals.max()))
    pad = 0.1 * (ymax - ymin)
    y_lo, y_hi = (ymin - pad, ymax + pad * 0.6)
    draw_prjeb55916(ax_left, activity, m55)
    draw_prjeb37009(ax_right, activity, m37)
    for ax in (ax_left, ax_right):
        ax.set_ylim(y_lo, y_hi)
    plt.setp(ax_right.get_yticklabels(), visible=False)
    ax_right.tick_params(axis="y", length=0)
    ax_right.spines["left"].set_visible(False)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=COL_CTRL,
            markersize=8,
            markeredgecolor="white",
            markeredgewidth=0.5,
            label="8 mM Q (control)",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=COL_PERT,
            markersize=8,
            markeredgecolor="white",
            markeredgewidth=0.5,
            label="0 mM Q (Gln deprivation)",
        ),
        Line2D(
            [0, 1],
            [0, 0],
            color="#444444",
            linewidth=2.2,
            linestyle="-",
            label="group mean",
        ),
    ]
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=5,
        frameon=False,
        fontsize=TICK_FONT,
        bbox_to_anchor=(0.5, -0.1),
    )
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    out_png = OUT_DIR / "figure4EF_projected_ATF4.png"
    out_pdf = OUT_DIR / "figure4EF_projected_ATF4.pdf"
    fig.savefig(out_png, dpi=DPI, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] saved -> {out_png}")
    print(f"[OK] saved -> {out_pdf}")


def load_inputs():
    raw = pd.read_excel(
        DATA / "validation.xlsx",
        sheet_name="Supplementary Table 15",
        header=None,
        nrows=107,
    )
    ids = raw.iloc[1, 1:].astype(str).str.strip()
    values = raw.iloc[2:, 1:].to_numpy(float)
    component_ids = raw.iloc[2:, 0].astype(int)
    activities = pd.DataFrame(values, index=component_ids, columns=ids)
    metadata = pd.read_excel(
        DATA / "validation.xlsx", sheet_name="Supplementary Table 14"
    )
    metadata["SampleID"] = metadata.SampleID.astype(str).str.strip()
    assert activities.columns.is_unique and metadata.SampleID.is_unique
    assert set(activities.columns) == set(metadata.SampleID)
    return activities, metadata


def project_meta(activities, metadata, project):
    selected = metadata[metadata.Project.eq(project)].copy()
    assert set(selected.SampleID) <= set(activities.columns)
    selected["col"] = selected.SampleID
    return selected


if __name__ == "__main__":
    activities, metadata = load_inputs()
    plot(activities, metadata)
