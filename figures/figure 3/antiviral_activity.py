"""Discovery stimulus groups; no mechanistic schematic.

Adapted from figure3D_antiviral_activity.py. Original inputs remain unchanged.
"""

import matplotlib

matplotlib.use("Agg")
from figure3_paths import DATA, RESULTS
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

PROJECT_DIR = DATA
XL_ACTIVITY = DATA / "discovery.xlsx"
XL_META = DATA / "metadata.xlsx"
OUT_DIR = RESULTS / "process_response"
OUT_DIR.mkdir(parents=True, exist_ok=True)
TARGET_IM = 87
TARGET_IM_NAME = "antiviral / dsRNA response"
ANTIVIRAL_PROJECT = "2015_12_FDA_virus"
VIRUS_PREFIXES = ["reo1", "reo3", "vsv1", "vsv2", "emcv1", "emcv2", "mvm1", "mvm2"]
VIRUS_COLORS = {"Reo": "#d62728", "VSV": "#1f77b4", "EMCV": "#ff7f0e", "MVM": "#9467bd"}
VIRUS_LABELS = {
    "Reo": "Reo (dsRNA virus)",
    "VSV": "VSV (ssRNA-)",
    "EMCV": "EMCV (ssRNA+)",
    "MVM": "MVM (ssDNA)",
}
VIRUS_ORDER = ["Reo", "VSV", "EMCV", "MVM"]
STIM_ORDER = ["mock", "virus_only", "ifn", "ifn_virus", "dsRNA", "dsRNA_virus"]
STIM_LABELS = [
    "\nMock\n(media)",
    "\nVirus\nalone",
    "\nIFN\nalone",
    "\nIFN +\nVirus",
    "\nPoly I:C\nalone",
    "\nPoly I:C\n+ Virus",
]
COL_LABELS = ["\nControl\n(other)"] + STIM_LABELS
COL_X = list(range(len(COL_LABELS)))
BASELINE_COLOR = "#9aa0a6"
BOX_FACE = "#e6e7ea"
BOX_EDGE = "#404040"
FIG_W, FIG_H = (4.6, 3.2)
DPI = 300
JITTER_WIDTH = 0.13
DOT_SIZE = 42
DOT_ALPHA = 0.92
DOT_EDGE = "white"
DOT_EDGE_W = 0.7
BASELINE_DOT_SIZE = 6
BASELINE_DOT_ALPHA = 0.18
BOX_WIDTH = 0.55
MEAN_BAR_HALF = 0.27
MEAN_BAR_LW = 2.2
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
    result.insert(0, "Condition", frame["Condition "].astype(str))
    assert result.SampleID.is_unique
    return result


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


def classify_stimulus(cond):
    sl = str(cond).lower()
    for prefix in VIRUS_PREFIXES:
        if sl.startswith(prefix):
            sl = sl[len(prefix) :]
            break
    has_virus = sl.startswith("v")
    if has_virus:
        sl = sl[1:]
    if "poly" in sl:
        return "dsRNA_virus" if has_virus else "dsRNA"
    if "ifn" in sl:
        return "ifn_virus" if has_virus else "ifn"
    if sl in ("m", "media", ""):
        return "virus_only" if has_virus else "mock"
    return "other"


def classify_virus_family(cond):
    sl = str(cond).lower()
    if sl.startswith("reo"):
        return "Reo"
    if sl.startswith("vsv"):
        return "VSV"
    if sl.startswith("emcv"):
        return "EMCV"
    if sl.startswith("mvm"):
        return "MVM"
    return None


def build_groups(A, im_col):
    df = A[["Project", "SampleID", "Condition", im_col]].copy()
    df = df.rename(columns={im_col: "activity"}).dropna(subset=["activity"])
    is_anti = df["Project"].astype(str) == ANTIVIRAL_PROJECT
    other = df[~is_anti].copy()
    fda = df[is_anti].copy()
    fda["arm"] = fda["Condition"].apply(classify_stimulus)
    fda["virus"] = fda["Condition"].apply(classify_virus_family)
    groups = {arm: fda[fda["arm"] == arm].copy() for arm in STIM_ORDER}
    groups["other"] = other
    return groups


def _style_box(bp, edgecolor=BOX_EDGE, facecolor=BOX_FACE):
    for box in bp["boxes"]:
        box.set(facecolor=facecolor, edgecolor=edgecolor, linewidth=1.1)
    for whisker in bp["whiskers"]:
        whisker.set(color=edgecolor, linewidth=1.0)
    for cap in bp["caps"]:
        cap.set(color=edgecolor, linewidth=1.0)
    for median in bp["medians"]:
        median.set(color="black", linewidth=1.4)
    for flier in bp.get("fliers", []):
        flier.set(marker="")


def plot_panel(groups, out_stem):
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    g0 = groups["other"]["activity"].to_numpy()
    bp = ax.boxplot(
        [g0],
        positions=[COL_X[0]],
        widths=BOX_WIDTH,
        patch_artist=True,
        showfliers=False,
        manage_ticks=False,
    )
    _style_box(bp, edgecolor=BOX_EDGE, facecolor=BOX_FACE)
    xs0 = COL_X[0] + RNG.uniform(-JITTER_WIDTH, JITTER_WIDTH, size=len(g0))
    ax.scatter(
        xs0,
        g0,
        s=BASELINE_DOT_SIZE,
        color=BASELINE_COLOR,
        alpha=BASELINE_DOT_ALPHA,
        edgecolor="none",
        zorder=2,
    )
    stim_means = {}
    for i, stim in enumerate(STIM_ORDER):
        x_pos = COL_X[i + 1]
        col_data = groups[stim]
        stim_means[stim] = (
            float(col_data["activity"].mean()) if len(col_data) else float("nan")
        )
        for fam in VIRUS_ORDER:
            sub = col_data[col_data["virus"] == fam]["activity"].to_numpy()
            if len(sub) == 0:
                continue
            xs = x_pos + RNG.uniform(-JITTER_WIDTH, JITTER_WIDTH, size=len(sub))
            ax.scatter(
                xs,
                sub,
                s=DOT_SIZE,
                color=VIRUS_COLORS[fam],
                alpha=DOT_ALPHA,
                edgecolor=DOT_EDGE,
                linewidth=DOT_EDGE_W,
                zorder=3,
            )
        m = stim_means[stim]
        if not np.isnan(m):
            ax.plot(
                [x_pos - MEAN_BAR_HALF, x_pos + MEAN_BAR_HALF],
                [m, m],
                color="black",
                linewidth=MEAN_BAR_LW,
                solid_capstyle="butt",
                zorder=4,
            )
    ax.set_xticks(COL_X)
    ax.set_xticklabels(COL_LABELS, fontsize=8.0)
    ax.set_xlim(-0.55, COL_X[-1] + 0.55)
    ax.set_ylabel(f"iM{TARGET_IM} activity", fontsize=10)
    ax.set_title(
        f"Interferon-Stimulated Genes (ISG)-1 (iM{TARGET_IM})", fontsize=10, pad=8
    )
    ax.grid(axis="y", linestyle=":", alpha=0.35)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="both", labelsize=8)
    y_low = ax.get_ylim()[0]
    n_per_col = [len(g0)] + [len(groups[s]) for s in STIM_ORDER]
    for x, n in zip(COL_X, n_per_col):
        ax.text(
            x,
            y_low - 0.5,
            f"n = {n}",
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
            markerfacecolor=VIRUS_COLORS[f],
            markeredgecolor="white",
            markersize=7,
            label=VIRUS_LABELS[f],
        )
        for f in VIRUS_ORDER
    ]
    handles.append(Line2D([0], [0], color="black", linewidth=MEAN_BAR_LW, label="mean"))
    fig.tight_layout()
    fig.savefig(out_stem.with_suffix(".png"), dpi=DPI, bbox_inches="tight")
    fig.savefig(out_stem.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    return stim_means


def main():
    print(f"[load] activity:  {XL_ACTIVITY}")
    A = load_activity_matrix(XL_ACTIVITY)
    print(f"[load] A: {A.shape}")
    im_col = f"iM{TARGET_IM}"
    if im_col not in A.columns:
        raise SystemExit(f"{im_col} not in A matrix")
    groups = build_groups(A, im_col)
    print("\n[groups]")
    other = groups["other"]
    print(
        f"  other        : n={len(other):3d}  median={float(np.median(other['activity'])):6.2f}  IQR=[{float(np.percentile(other['activity'], 25)):.2f}, {float(np.percentile(other['activity'], 75)):.2f}]  95th %ile={float(np.percentile(other['activity'], 95)):.2f}"
    )
    for k in STIM_ORDER:
        v = groups[k]
        mean = v["activity"].mean() if len(v) else float("nan")
        print(f"  {k:13s}: n={len(v):2d}  mean={mean:6.2f}")
    out_stem = OUT_DIR / "figure3D_antiviral_activity"
    stim_means = plot_panel(groups, out_stem)
    print(f"\n[save] {out_stem}.png")
    print(f"[save] {out_stem}.pdf")
    mock = groups["mock"]["activity"].to_numpy()
    other_arr = groups["other"]["activity"].to_numpy()
    print(f"\n[Hedges g vs mock (n=8)]")
    for stim in STIM_ORDER:
        if stim == "mock":
            continue
        v = groups[stim]["activity"].to_numpy()
        print(f"  {stim:13s}  g = {hedges_g(mock, v):5.2f}")
    print(f"\n[Hedges g vs other-projects baseline (n={len(other_arr)})]")
    for stim in STIM_ORDER:
        v = groups[stim]["activity"].to_numpy()
        print(f"  {stim:13s}  g = {hedges_g(other_arr, v):5.2f}")


if __name__ == "__main__":
    main()
