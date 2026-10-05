"""Discovery TDS activity and Hedges g.

Adapted from figure3A_temperature_activity.py. Original inputs remain unchanged.
"""

import matplotlib

matplotlib.use("Agg")
from figure3_paths import DATA, RESULTS
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PROJECT_DIR = DATA
XL_ACTIVITY = DATA / "discovery.xlsx"
XL_META = DATA / "metadata.xlsx"
OUT_DIR = RESULTS / "process_response"
OUT_DIR.mkdir(parents=True, exist_ok=True)
TARGET_IM = 39
TARGET_IM_NAME = "Temperature downshift"
TS_PROJECTS = ["PRJNA593052", "GC"]
PROJECT_COLORS = {"PRJNA593052": "#3b82bd", "GC": "#f5a623"}
PROJECT_LABELS = {
    "PRJNA593052": "PRJNA593052 (CHO-K1, 72 h)",
    "GC": "GC (CHO-DG44 / -K1, Day 9)",
}
COL_KEYS = ["other_no", "matched_no", "ts_yes"]
COL_LABELS = [
    "\nControl\n(other projects)",
    "\nControl\n(matched)",
    "\nHypothermia\n(32 C)",
]
COL_X = [0.0, 1.0, 2.0]
BASELINE_COLOR = "#9aa0a6"
FIG_W, FIG_H = (4.6, 3.2)
DPI = 300
JITTER_WIDTH = 0.1
DOT_SIZE = 36
DOT_ALPHA = 0.85
DOT_EDGE = "white"
DOT_EDGE_W = 0.6
MEAN_MARKER_SIZE = 110
MEAN_LINE_W = 2.0
BASELINE_DOT_SIZE = 7
BASELINE_DOT_ALPHA = 0.25
VIOLIN_WIDTH = 0.55
VIOLIN_ALPHA = 0.35
RNG = np.random.default_rng(20260426)


def load_activity_matrix(xl_path):
    frame = pd.read_excel(xl_path, sheet_name="A_reorganized")
    columns = [c for c in frame if isinstance(c, (int, np.integer))]
    result = (
        frame[columns]
        .apply(pd.to_numeric, errors="raise")
        .rename(columns=lambda c: f"iM{c}")
    )
    result.insert(0, "SampleID", frame["Unnamed: 1"].astype(str))
    result.insert(0, "Project", frame["Unnamed: 0"].astype(str))
    assert result.SampleID.is_unique
    return result


def load_metadata(xl_path):
    md = pd.read_excel(xl_path)
    md["Temperature shift"] = md["Temperature shift"].astype(str).str.strip()
    return md


def hedges_g(x, y):
    nx, ny = (len(x), len(y))
    if nx < 2 or ny < 2:
        return float("nan")
    sx2, sy2 = (np.var(x, ddof=1), np.var(y, ddof=1))
    sp = np.sqrt(((nx - 1) * sx2 + (ny - 1) * sy2) / (nx + ny - 2))
    if sp == 0:
        return float("nan")
    d = (np.mean(y) - np.mean(x)) / sp
    j = 1 - 3 / (4 * (nx + ny) - 9)
    return float(d * j)


def build_groups(A, md, im_col, ts_projects):
    """Return dict with three groups + per-project subgroups for Cols 2 and 3."""
    df = (
        md[["SampleID", "Project", "Temperature shift"]]
        .merge(A[["SampleID", im_col]], on="SampleID", how="left")
        .dropna(subset=[im_col])
    )
    df = df.rename(columns={im_col: "activity", "Temperature shift": "TempShift"})
    other_no = df[~df["Project"].isin(ts_projects) & (df["TempShift"] == "NO")]
    matched_no = df[df["Project"].isin(ts_projects) & (df["TempShift"] == "NO")]
    ts_yes = df[df["Project"].isin(ts_projects) & (df["TempShift"] == "YES")]
    return {"other_no": other_no, "matched_no": matched_no, "ts_yes": ts_yes}


def plot_panel(groups, out_stem):
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    g1 = groups["other_no"]["activity"].to_numpy()
    g2_all = groups["matched_no"]["activity"].to_numpy()
    g3_all = groups["ts_yes"]["activity"].to_numpy()
    box = ax.boxplot(
        [g1, g2_all, g3_all],
        positions=COL_X,
        widths=VIOLIN_WIDTH,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=1.6),
        boxprops=dict(linewidth=1.0, color="0.35"),
        whiskerprops=dict(linewidth=1.0, color="0.35"),
        capprops=dict(linewidth=1.0, color="0.35"),
    )
    for patch in box["boxes"]:
        patch.set_facecolor("white")
        patch.set_alpha(0.95)
    xs = COL_X[0] + RNG.uniform(-JITTER_WIDTH, JITTER_WIDTH, size=len(g1))
    ax.scatter(
        xs,
        g1,
        s=BASELINE_DOT_SIZE,
        color=BASELINE_COLOR,
        alpha=BASELINE_DOT_ALPHA,
        edgecolor="none",
        zorder=2,
    )
    summaries = {}
    for proj, color in PROJECT_COLORS.items():
        g2 = groups["matched_no"]
        g2 = g2[g2["Project"] == proj]["activity"].to_numpy()
        g3 = groups["ts_yes"]
        g3 = g3[g3["Project"] == proj]["activity"].to_numpy()
        if len(g2):
            xs = COL_X[1] + RNG.uniform(-JITTER_WIDTH, JITTER_WIDTH, size=len(g2))
            ax.scatter(
                xs,
                g2,
                s=DOT_SIZE,
                color=color,
                alpha=DOT_ALPHA,
                edgecolor=DOT_EDGE,
                linewidth=DOT_EDGE_W,
                zorder=3,
            )
        if len(g3):
            xs = COL_X[2] + RNG.uniform(-JITTER_WIDTH, JITTER_WIDTH, size=len(g3))
            ax.scatter(
                xs,
                g3,
                s=DOT_SIZE,
                color=color,
                alpha=DOT_ALPHA,
                edgecolor=DOT_EDGE,
                linewidth=DOT_EDGE_W,
                zorder=3,
            )
        delta = g3.mean() - g2.mean() if len(g2) and len(g3) else float("nan")
        g = hedges_g(g2, g3)
        summaries[proj] = {
            "n_matched_no": int(len(g2)),
            "n_ts_yes": int(len(g3)),
            "mean_matched_no": float(g2.mean()) if len(g2) else float("nan"),
            "mean_ts_yes": float(g3.mean()) if len(g3) else float("nan"),
            "delta": float(delta),
            "hedges_g": float(g),
        }
    med1 = float(np.median(g1))
    q25, q75 = np.percentile(g1, [25, 75])
    ax.set_xticks(COL_X)
    ax.set_xticklabels(COL_LABELS, fontsize=8.5)
    ax.set_xlim(-0.55, 2.55)
    ax.set_ylabel(f"iM{TARGET_IM} activity", fontsize=10)
    ax.set_title(f"Temperature Downshift (iM{TARGET_IM})", fontsize=10, pad=8)
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="both", labelsize=8)
    n_other = len(groups["other_no"])
    n_matched = len(groups["matched_no"])
    n_ts = len(groups["ts_yes"])
    n_text = [f"n = {n_other}", f"n = {n_matched}", f"n = {n_ts}"]
    for x, t in zip(COL_X, n_text):
        ax.text(
            x,
            ax.get_ylim()[0] - 0.5,
            t,
            ha="center",
            va="top",
            fontsize=7.5,
            color="0.35",
        )
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=color,
            markeredgecolor="white",
            markersize=7,
            label=PROJECT_LABELS[p],
        )
        for p, color in PROJECT_COLORS.items()
        if p in summaries
    ]
    if handles:
        ax.legend(
            handles=handles,
            loc="upper left",
            fontsize=7,
            frameon=True,
            framealpha=0.95,
            borderpad=0.4,
            handletextpad=0.4,
        )
    fig.tight_layout()
    fig.savefig(out_stem.with_suffix(".png"), dpi=DPI, bbox_inches="tight")
    fig.savefig(out_stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return (
        summaries,
        {
            "baseline_median": med1,
            "baseline_q25": q25,
            "baseline_q75": q75,
            "n_other_no": n_other,
            "n_matched_no": n_matched,
            "n_ts_yes": n_ts,
        },
    )


def main():
    print(f"[load] activity:  {XL_ACTIVITY}")
    print(f"[load] metadata:  {XL_META}")
    A = load_activity_matrix(XL_ACTIVITY)
    md = load_metadata(XL_META)
    print(f"[load] A: {A.shape}, md: {md.shape}")
    im_col = f"iM{TARGET_IM}"
    if im_col not in A.columns:
        raise SystemExit(f"{im_col} not in A matrix")
    groups = build_groups(A, md, im_col, TS_PROJECTS)
    print("\n[groups]")
    for k, v in groups.items():
        print(f"  {k}: n={len(v)}")
    out_stem = OUT_DIR / "figure3A_temperature_activity"
    summaries, baseline = plot_panel(groups, out_stem)
    print(f"\n[save] {out_stem}.png")
    print(f"[save] {out_stem}.pdf")
    print("\n[baseline]")
    print(
        f"  Col 1 median = {baseline['baseline_median']:.2f}  IQR = [{baseline['baseline_q25']:.2f}, {baseline['baseline_q75']:.2f}]  n = {baseline['n_other_no']}"
    )
    print("\n[summary]")
    for proj, s in summaries.items():
        print(
            f"  {proj}: matched-NO n={s['n_matched_no']} (mean {s['mean_matched_no']:.2f}) -> TS n={s['n_ts_yes']} (mean {s['mean_ts_yes']:.2f})  D=+{s['delta']:.2f}  g={s['hedges_g']:.2f}"
        )


if __name__ == "__main__":
    main()
