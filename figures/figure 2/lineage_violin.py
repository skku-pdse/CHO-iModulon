"""Original manuscript violin design.

Adapted from figure4A_violin_plot.py. Original inputs remain unchanged.
"""

import matplotlib

matplotlib.use("Agg")
from figure2_paths import DATA, RESULTS
import re
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

try:
    import scikit_posthocs as sp

    HAS_POSTHOCS = True
except ImportError:
    raise ImportError(
        "Install scikit-posthocs from requirements.txt to reproduce Dunn tests"
    )
IMODULON_FILE = DATA / "discovery.xlsx"
METADATA_FILE = DATA / "metadata.xlsx"
OUTDIR = RESULTS / "lineage"
OUTDIR.mkdir(exist_ok=True)


def normalize_sample_name(x: str) -> str:
    if pd.isna(x):
        return x
    x = str(x).strip()
    x = re.sub("\\s+", "", x)
    return x


def load_a_matrix(imodulon_file: str, sheet_name: str = "A") -> pd.DataFrame:
    """
    A sheet
    row: iModulon
    col: sample
    """
    A = pd.read_excel(imodulon_file, sheet_name=sheet_name, index_col=0, header=1)
    A.columns = [normalize_sample_name(c) for c in A.columns]
    A.index = [str(i) for i in A.index]
    return A


def load_metadata(metadata_file: str, sheet_name: str = "Main_v3") -> pd.DataFrame:
    meta = pd.read_excel(metadata_file, sheet_name=sheet_name)
    if "SampleID" not in meta.columns:
        raise ValueError(f"'SampleID' column not found. columns={list(meta.columns)}")
    meta = meta.copy()
    meta["Sample"] = meta["SampleID"].map(normalize_sample_name)
    # Current activity and metadata inputs already share the corrected D3 identifier.
    if meta["Sample"].duplicated().any():
        counts = {}
        new_names = []
        for s in meta["Sample"]:
            counts[s] = counts.get(s, 0) + 1
            if counts[s] == 1:
                new_names.append(s)
            else:
                new_names.append(f"{s}.{counts[s] - 1}")
        meta["Sample"] = new_names
    return meta


def merge_a_and_metadata(A: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    A_long = (
        A.T.reset_index()
        .rename(columns={"index": "Sample"})
        .melt(id_vars="Sample", var_name="iModulon", value_name="Activity")
    )
    assert meta.Sample.is_unique and A.columns.is_unique
    assert set(A.columns) == set(meta.Sample), "Unmatched sample identifiers"
    merged = A_long.merge(meta, on="Sample", how="inner", validate="many_to_one")
    return merged


def p_to_stars(p):
    if pd.isna(p):
        return "ns"
    if p < 0.0001:
        return "****"
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    return "ns"


def get_group_order(df, metadata_col, group_order=None):
    if group_order is not None:
        return [g for g in group_order if g in df[metadata_col].astype(str).unique()]
    return df[metadata_col].astype(str).value_counts().index.tolist()


def build_xticklabels_with_n(df, metadata_col, order):
    counts = df.groupby(metadata_col)["Sample"].nunique()
    labels = []
    for g in order:
        n = int(counts.get(g, 0))
        labels.append(f"{g}\n(n={n})")
    return labels


def compute_pairwise_dunn(sub_df, metadata_col, order):
    """
    Dunn's post hoc with BH correction
    returns dataframe: group1, group2, qval
    """
    if not HAS_POSTHOCS:
        return pd.DataFrame(columns=["group1", "group2", "qval"])
    tmp = sub_df[[metadata_col, "Activity"]].copy()
    tmp[metadata_col] = tmp[metadata_col].astype(str)
    dunn = sp.posthoc_dunn(
        tmp, val_col="Activity", group_col=metadata_col, p_adjust="fdr_bh"
    )
    rows = []
    for g1, g2 in combinations(order, 2):
        if g1 in dunn.index and g2 in dunn.columns:
            rows.append({"group1": g1, "group2": g2, "qval": float(dunn.loc[g1, g2])})
    return pd.DataFrame(rows)


def add_sig_brackets(ax, pairs_df, order, y_start, y_step, linewidth=1.2, color="0.35"):
    pos_map = {g: i for i, g in enumerate(order)}
    current_y = y_start
    for _, row in pairs_df.iterrows():
        g1, g2, qv = (row["group1"], row["group2"], row["qval"])
        if g1 not in pos_map or g2 not in pos_map:
            continue
        x1, x2 = (pos_map[g1], pos_map[g2])
        if x1 > x2:
            x1, x2 = (x2, x1)
        ax.plot(
            [x1, x1, x2, x2],
            [
                current_y,
                current_y + y_step * 0.25,
                current_y + y_step * 0.25,
                current_y,
            ],
            lw=linewidth,
            c=color,
            clip_on=False,
        )
        ax.text(
            (x1 + x2) / 2,
            current_y + y_step * 0.3,
            p_to_stars(qv),
            ha="center",
            va="bottom",
            fontsize=11,
            fontweight="bold",
            color="black",
        )
        current_y += y_step


def plot_selected_imodulons_by_metadata(
    merged_df: pd.DataFrame,
    metadata_col: str,
    imodulons: list,
    group_order: list = None,
    title_map: dict = None,
    panel_colors: dict = None,
    panel_colors_violin: dict = None,
    outdir: Path = OUTDIR,
    save_name: str = None,
    min_group_size: int = 2,
    add_stats: bool = True,
    pairs_to_annotate: dict = None,
    figsize_per_panel=(4.2, 5.6),
    sharey=True,
):
    """Plot selected activities by metadata group with optional post-hoc tests."""
    if metadata_col not in merged_df.columns:
        raise ValueError(f"metadata column '{metadata_col}' not found.")
    df = merged_df.copy()
    df = df[df[metadata_col].notna()].copy()
    df[metadata_col] = df[metadata_col].astype(str)
    valid_groups = (
        df.groupby(metadata_col)["Sample"]
        .nunique()
        .loc[lambda x: x >= min_group_size]
        .index.tolist()
    )
    df = df[df[metadata_col].isin(valid_groups)].copy()
    if df.empty:
        raise ValueError(
            f"No usable data after filtering metadata_col='{metadata_col}'."
        )
    order = get_group_order(df, metadata_col, group_order)
    n = len(imodulons)
    ncols = min(4, n)
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows),
        sharey=sharey,
        constrained_layout=True,
    )
    if n == 1:
        axes = np.array([axes])
    axes = np.array(axes).reshape(-1)
    if title_map is None:
        title_map = {}
    if panel_colors is None:
        panel_colors = {}
    sns.set_theme(style="white", context="talk")
    for ax, imod in zip(axes, imodulons):
        imod = str(imod)
        sub = df[df["iModulon"].astype(str) == imod].copy()
        if sub.empty:
            ax.set_visible(False)
            continue
        color_violin = panel_colors_violin.get(imod, "#4c78a8")
        color = panel_colors.get(imod, "#4c78a8")
        sns.violinplot(
            data=sub,
            x=metadata_col,
            y="Activity",
            order=order,
            inner=None,
            cut=0,
            linewidth=1.4,
            color=color_violin,
            saturation=1,
            ax=ax,
        )
        sns.stripplot(
            data=sub,
            x=metadata_col,
            y="Activity",
            order=order,
            color=color,
            alpha=0.8,
            size=2.7,
            jitter=0.18,
            linewidth=0.5,
            ax=ax,
        )
        panel_title = title_map.get(imod, f"iM{imod}")
        ax.set_title(panel_title, fontsize=15, fontweight="bold", color=color, pad=16)
        ax.set_xlabel("")
        if ax == axes[0]:
            ax.set_ylabel("Activity", fontsize=18, fontweight="bold")
        else:
            ax.set_ylabel("")
        xticklabels = build_xticklabels_with_n(sub, metadata_col, order)
        ax.set_xticklabels(xticklabels, rotation=0, fontsize=10)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.axhline(0, color="0.75", linewidth=1, zorder=0)
        if add_stats:
            groups = [
                sub.loc[sub[metadata_col] == g, "Activity"].dropna().values
                for g in order
            ]
            groups = [g for g in groups if len(g) > 0]
            if len(groups) >= 2:
                try:
                    _ = kruskal(*groups)
                except Exception:
                    pass
            y_min = sub["Activity"].min()
            y_max = sub["Activity"].max()
            y_range = max(y_max - y_min, 1.0)
            y_start = y_max + y_range * 0.1
            y_step = y_range * 0.12
            pair_df = compute_pairwise_dunn(sub, metadata_col, order)
            if not pair_df.empty:
                if pairs_to_annotate is not None and imod in pairs_to_annotate:
                    wanted = set((tuple(x) for x in pairs_to_annotate[imod]))
                    pair_df = pair_df[
                        pair_df.apply(
                            lambda r: (r["group1"], r["group2"]) in wanted
                            or (r["group2"], r["group1"]) in wanted,
                            axis=1,
                        )
                    ].copy()
                else:
                    pair_df = (
                        pair_df[pair_df["qval"] < 0.05]
                        .sort_values("qval", ascending=True)
                        .head(4)
                    )
                if not pair_df.empty:
                    add_sig_brackets(ax, pair_df, order, y_start=y_start, y_step=y_step)
                    ax.set_ylim(
                        y_min - y_range * 0.08, y_start + y_step * (len(pair_df) + 0.7)
                    )
                else:
                    ax.set_ylim(y_min - y_range * 0.08, y_max + y_range * 0.18)
        ax.tick_params(axis="y", labelsize=11)
    for ax in axes[len(imodulons) :]:
        ax.set_visible(False)
    if save_name is None:
        save_name = f"{metadata_col.replace('/', '_').replace(' ', '_')}_{'_'.join(map(str, imodulons))}"
    out_png = outdir / f"{save_name}.png"
    out_pdf = outdir / f"{save_name}.pdf"
    fig.savefig(out_png, dpi=400, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close("all")
    print(f"[SAVED] {out_png}")
    print(f"[SAVED] {out_pdf}")


if __name__ == "__main__":
    A = load_a_matrix(IMODULON_FILE, sheet_name="A")
    meta = load_metadata(METADATA_FILE, sheet_name="Main_v3")
    merged = merge_a_and_metadata(A, meta)
    plot_selected_imodulons_by_metadata(
        merged_df=merged,
        metadata_col="Cell line",
        imodulons=["52", "65"],
        group_order=["CHO-K1", "CHO-DG44", "CHO-S", "CHO-ZeLa"],
        title_map={"52": "ECM-1 (iM52)", "65": "ECM-2 (iM65)"},
        panel_colors_violin={"65": "#8FB6E0", "52": "#E89A8D"},
        panel_colors={"65": "#2f69b0", "52": "#c0392b"},
        save_name="figure2AB_lineage_activity",
        min_group_size=2,
        add_stats=True,
        pairs_to_annotate=None,
    )
