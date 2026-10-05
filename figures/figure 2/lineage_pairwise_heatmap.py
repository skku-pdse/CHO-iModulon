"""Selected pairwise display; original filename retained in outputs.

Adapted from figure4B_lineage_pairwise_heatmap_1x3.py. Original inputs remain unchanged.
"""

import matplotlib

matplotlib.use("Agg")
from figure2_paths import DATA, RESULTS
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cm
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors

PANEL_DIR = RESULTS / "lineage"
IMODULON_FILE = DATA / "discovery.xlsx"
OUTDIR = RESULTS / "lineage"
OUTDIR.mkdir(parents=True, exist_ok=True)
GROUP_ORDER = ["CHO-K1", "CHO-DG44", "CHO-S", "CHO-ZeLa"]
IMODULON_LIST = ["52", "65"]
PANEL_SIZE = 3.4
VMAX_D = 3.0
base_1 = plt.cm.RdBu_r
CMAP_D = mcolors.LinearSegmentedColormap.from_list(
    "trunc_oranges", base_1(np.linspace(0.05, 0.95, 256))
)
base_2 = plt.cm.Oranges
CMAP_AUC = mcolors.LinearSegmentedColormap.from_list(
    "trunc_oranges", base_2(np.linspace(0.0, 0.75, 256))
)
ANNOT_FONTSIZE = 10
TICK_FONTSIZE = 10
TITLE_FONTSIZE = 11
SUPTITLE_FONTSIZE = 12
DIAG_HATCH = ""
DIAG_FACE = "white"
DIAG_EDGE = "0.7"


def load_inputs():
    pairwise = pd.read_csv(PANEL_DIR / "lineage_pairwise_effect.csv")
    summary = pd.read_csv(PANEL_DIR / "lineage_effect_summary.csv")
    pairwise["iModulon"] = pairwise["iModulon"].astype(str)
    summary["iModulon"] = summary["iModulon"].astype(str)
    try:
        names = pd.read_excel(IMODULON_FILE, sheet_name="iModulon Summary")
        names["iModulon"] = names["iModulon"].astype(str)
        name_map = dict(zip(names["iModulon"], names["Name"]))
    except Exception as e:
        print(f"[WARN] Could not load iModulon names: {e}")
        name_map = {}
    return (pairwise, summary, name_map)


def build_matrices(pairwise_df, imod, group_order=GROUP_ORDER):
    """
    Returns (d_mat, auc_mat) as 4x4 DataFrames indexed/columned by group_order.

    Convention:
      d_mat[row, col] = mean(row) - mean(col), standardized (Cohen's d)
      auc_mat[row, col] = |AUROC(row vs col) - 0.5| * 2

    d_mat is anti-symmetric, auc_mat is symmetric.
    """
    sub = pairwise_df[pairwise_df["iModulon"] == imod]
    n = len(group_order)
    d = pd.DataFrame(np.nan, index=group_order, columns=group_order, dtype=float)
    a = pd.DataFrame(np.nan, index=group_order, columns=group_order, dtype=float)
    for _, r in sub.iterrows():
        ga, gb = (r["group_A"], r["group_B"])
        if ga not in group_order or gb not in group_order:
            continue
        d.loc[ga, gb] = r["cohens_d"]
        d.loc[gb, ga] = -r["cohens_d"]
        auc_disc = r["auroc_disc"] * 2 if not np.isnan(r["auroc_disc"]) else np.nan
        a.loc[ga, gb] = auc_disc
        a.loc[gb, ga] = auc_disc
    return (d, a)


def plot_one_panel(ax, d_mat, auc_mat, title, vmax_d=VMAX_D):
    """
    Lower triangle (i > j): signed Cohen's d, diverging cmap.
    Upper triangle (i < j): AUROC discrimination, sequential cmap.
    Diagonal: hatched.
    """
    n = d_mat.shape[0]
    groups = d_mat.index.tolist()
    d_lower = np.full((n, n), np.nan)
    auc_upper = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            if i > j:
                d_lower[i, j] = d_mat.iloc[i, j]
            elif i < j:
                auc_upper[i, j] = auc_mat.iloc[i, j]
    norm_d = Normalize(vmin=-vmax_d, vmax=vmax_d, clip=True)
    norm_a = Normalize(vmin=0, vmax=1, clip=True)
    cmap_d = plt.get_cmap(CMAP_D)
    cmap_a = plt.get_cmap(CMAP_AUC)
    for i in range(n):
        for j in range(n):
            x = j
            y = n - 1 - i
            if i > j:
                v = d_lower[i, j]
                if np.isnan(v):
                    color = (1, 1, 1, 0)
                    fmt = ""
                    txt_color = "black"
                else:
                    color = cmap_d(norm_d(v))
                    txt_color = "white" if abs(v) > vmax_d * 0.6 else "black"
                    if abs(v) >= vmax_d * 1.5:
                        fmt = f"{v:.0f}"
                    else:
                        fmt = f"{v:.2f}"
                ax.add_patch(
                    Rectangle(
                        (x, y), 1, 1, facecolor=color, edgecolor="0.7", linewidth=0.6
                    )
                )
                if fmt:
                    ax.text(
                        x + 0.5,
                        y + 0.5,
                        fmt,
                        ha="center",
                        va="center",
                        fontsize=ANNOT_FONTSIZE,
                        color=txt_color,
                    )
            elif i < j:
                v = auc_upper[i, j]
                if np.isnan(v):
                    color = (1, 1, 1, 0)
                    fmt = ""
                    txt_color = "black"
                else:
                    color = cmap_a(norm_a(v))
                    txt_color = "white" if v > 0.6 else "black"
                    fmt = f"{v:.2f}"
                ax.add_patch(
                    Rectangle(
                        (x, y), 1, 1, facecolor=color, edgecolor="0.7", linewidth=0.6
                    )
                )
                if fmt:
                    ax.text(
                        x + 0.5,
                        y + 0.5,
                        fmt,
                        ha="center",
                        va="center",
                        fontsize=ANNOT_FONTSIZE,
                        color=txt_color,
                    )
            else:
                ax.add_patch(
                    Rectangle(
                        (x, y),
                        1,
                        1,
                        facecolor=DIAG_FACE,
                        edgecolor=DIAG_EDGE,
                        linewidth=0.6,
                        hatch=DIAG_HATCH,
                    )
                )
    ax.set_xlim(0, n)
    ax.set_ylim(0, n)
    ax.set_xticks(np.arange(n) + 0.5)
    ax.set_yticks(np.arange(n) + 0.5)
    short = [g.replace("CHO-", "") for g in groups]
    ax.set_xticklabels(short, fontsize=TICK_FONTSIZE, rotation=0)
    ax.set_yticklabels(short[::-1], fontsize=TICK_FONTSIZE)
    ax.set_aspect("equal")
    ax.set_title(title, fontsize=TITLE_FONTSIZE, pad=8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(length=0)


def make_figure(
    pairwise_df, summary_df, name_map, imod_list, panel_size=PANEL_SIZE, vmax_d=VMAX_D
):
    n_panels = len(imod_list)
    n_cols = n_panels
    n_rows = 1
    fig_w = n_cols * panel_size + 0.6
    fig_h = panel_size + 1.6
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(2, n_cols, height_ratios=[1, 0.1], hspace=0.55, wspace=0.32)
    summary_idx = summary_df.set_index("iModulon")
    for k, imod in enumerate(imod_list):
        ax = fig.add_subplot(gs[0, k])
        d_mat, auc_mat = build_matrices(pairwise_df, imod)
        nm = name_map.get(imod, "")
        if imod in summary_idx.index:
            row = summary_idx.loc[imod]
            if "omega_sq" in row.index and (not pd.isna(row["omega_sq"])):
                eff = row["omega_sq"]
                eff_label = "ω²"
            else:
                eff = row["eta_sq"]
                eff_label = "η²"
            bim = bool(row.get("any_bimodal_flag", False))
            outlier = str(row.get("outlier_lineage", ""))
            direction = str(row.get("direction", ""))
            bim_mark = " *" if bim else ""
            outlier_short = (
                outlier.replace("CHO-", "") if outlier and outlier != "nan" else "—"
            )
            dir_short = direction if direction and direction != "nan" else ""
            sub = f"{eff_label}={eff:.2f} | {outlier_short} {dir_short}"
            title = f"{nm} (iM{imod})\n{sub}"
        else:
            title = f"{nm} (iM{imod})"
        plot_one_panel(ax, d_mat, auc_mat, title, vmax_d=vmax_d)
    if n_cols >= 2:
        cax_d = fig.add_subplot(gs[1, 0])
        cax_a = fig.add_subplot(gs[1, n_cols - 1])
    else:
        cax_d = fig.add_subplot(gs[1, 0])
        cax_a = fig.add_subplot(gs[1, 0])
    norm_d = Normalize(vmin=-vmax_d, vmax=vmax_d)
    norm_a = Normalize(vmin=0, vmax=1)
    cb_d = plt.colorbar(
        cm.ScalarMappable(norm=norm_d, cmap=CMAP_D),
        cax=cax_d,
        orientation="horizontal",
        extend="both",
    )
    cb_d.set_label(
        f"Lower △: Cohen's d (row − col), saturated at ±{vmax_d}", fontsize=9
    )
    cb_d.ax.tick_params(labelsize=8)
    cb_a = plt.colorbar(
        cm.ScalarMappable(norm=norm_a, cmap=CMAP_AUC),
        cax=cax_a,
        orientation="horizontal",
    )
    cb_a.set_label("Upper △: AUROC discrimination |AUC−0.5|×2", fontsize=9)
    cb_a.ax.tick_params(labelsize=8)
    return fig


if __name__ == "__main__":
    print("[INFO] Loading inputs ...")
    pairwise_df, summary_df, name_map = load_inputs()
    imod_list = [str(x) for x in IMODULON_LIST]
    missing = [im for im in imod_list if im not in summary_df["iModulon"].values]
    if missing:
        print(f"[WARN] Selected iModulons not found in summary: {missing}")
    print(f"[INFO] Selected iModulons: {imod_list}")
    print("[INFO] Building 1x3 figure ...")
    fig = make_figure(pairwise_df, summary_df, name_map, imod_list)
    out_png = OUTDIR / "figure2CD_lineage_pairwise.png"
    out_pdf = OUTDIR / "figure2CD_lineage_pairwise.pdf"
    fig.savefig(out_png, dpi=400, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    print(f"[SAVED] {out_png}")
    print(f"[SAVED] {out_pdf}")
