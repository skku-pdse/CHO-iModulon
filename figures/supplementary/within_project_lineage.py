"""Current multi-project discovery activity comparison.

Adapted from supplementary_figure2_multilineage.py. Original inputs remain unchanged.
"""

import matplotlib

matplotlib.use("Agg")
from supplementary_paths import DATA, RESULTS
from pathlib import Path
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

ROOT = DATA
ACTIVITY_XLSX = DATA / "discovery.xlsx"
METADATA_XLSX = DATA / "metadata.xlsx"
OUT = RESULTS / "lineage"
PROJECT_ORDER = [
    "GC",
    "PRJNA318886",
    "2014_11_CHOS_baseline",
    "2020_02_12_ZeN",
    "2016_05_HH_L",
    "2017_10_ZeLa_C13",
]
PROJECT_LABELS = {
    "GC": "GC",
    "PRJNA318886": "PRJNA318886",
    "2014_11_CHOS_baseline": "2014_11_CHOS\nbaseline",
    "2020_02_12_ZeN": "2020_02_12_ZeN",
    "2016_05_HH_L": "2016_05_HH_L",
    "2017_10_ZeLa_C13": "2017_10_ZeLa\nC13",
}
MODULES = [
    (52, "ECM-1 (discovery iM52)", "#C83B2D"),
    (65, "ECM-2 (discovery iM65)", "#2F6DB2"),
]


def load_data() -> pd.DataFrame:
    activity_raw = pd.read_excel(ACTIVITY_XLSX, sheet_name="A_reorganized")
    sample_col = "Unnamed: 1"
    module_cols = [
        col for col in activity_raw.columns if isinstance(col, (int, np.integer))
    ]
    activity = (
        activity_raw.dropna(subset=[sample_col])
        .drop_duplicates(subset=sample_col, keep="first")
        .assign(**{sample_col: lambda frame: frame[sample_col].astype(str)})
        .set_index(sample_col)[module_cols]
        .apply(pd.to_numeric, errors="coerce")
    )
    metadata = pd.read_excel(METADATA_XLSX, sheet_name="Main_v3")
    metadata.columns = metadata.columns.astype(str).str.strip()
    metadata["SampleID"] = metadata["SampleID"].astype(str)
    metadata = metadata.drop_duplicates("SampleID", keep="first").set_index("SampleID")
    common = activity.index.intersection(metadata.index)
    data = metadata.loc[common, ["Project", "Cell line"]].join(
        activity.loc[common, [52, 65]]
    )
    data = data[data["Project"].isin(PROJECT_ORDER)].copy()
    return data


def omega_squared(values: pd.Series, groups: pd.Series) -> float:
    frame = pd.DataFrame({"y": values, "group": groups}).dropna()
    n = len(frame)
    levels = frame["group"].unique()
    k = len(levels)
    if k < 2 or n <= k:
        return np.nan
    grand = frame["y"].mean()
    ss_between = sum(
        (
            len(group) * (group["y"].mean() - grand) ** 2
            for _, group in frame.groupby("group", observed=True)
        )
    )
    ss_total = ((frame["y"] - grand) ** 2).sum()
    ss_within = ss_total - ss_between
    ms_within = ss_within / (n - k)
    return max((ss_between - (k - 1) * ms_within) / (ss_total + ms_within), 0)


def draw() -> None:
    data = load_data()
    OUT.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="ticks", context="paper")
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 9,
            "axes.labelsize": 10,
            "axes.titlesize": 10,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(
        2,
        len(PROJECT_ORDER),
        figsize=(14.2, 6.4),
        sharey="row",
        gridspec_kw={"hspace": 0.3, "wspace": 0.14},
    )
    rng = np.random.default_rng(11)
    for row, (module, module_name, module_color) in enumerate(MODULES):
        for col, project in enumerate(PROJECT_ORDER):
            ax = axes[row, col]
            subset = data[data["Project"].eq(project)].dropna(subset=[module]).copy()
            lineages = sorted(subset["Cell line"].dropna().unique())
            sns.boxplot(
                data=subset,
                x="Cell line",
                y=module,
                order=lineages,
                width=0.52,
                showfliers=False,
                color="white",
                boxprops={"edgecolor": module_color, "linewidth": 1.25},
                medianprops={"color": module_color, "linewidth": 1.5},
                whiskerprops={"color": module_color, "linewidth": 1.0},
                capprops={"color": module_color, "linewidth": 1.0},
                ax=ax,
            )
            for xpos, lineage in enumerate(lineages):
                values = subset.loc[subset["Cell line"].eq(lineage), module].to_numpy()
                jitter = rng.uniform(-0.16, 0.16, size=len(values))
                ax.scatter(
                    xpos + jitter,
                    values,
                    s=13,
                    color=module_color,
                    alpha=0.62,
                    edgecolor="white",
                    linewidth=0.25,
                    zorder=3,
                )
            effect = omega_squared(subset[module], subset["Cell line"])
            ax.text(
                0.97,
                0.05,
                f"$\\omega^2$ = {effect:.2f}",
                transform=ax.transAxes,
                ha="right",
                va="bottom",
                fontsize=8,
                color="#333333",
            )
            counts = subset["Cell line"].value_counts()
            ax.set_xticks(range(len(lineages)))
            ax.set_xticklabels(
                [f"{lineage}\n(n={counts[lineage]})" for lineage in lineages],
                rotation=0,
                ha="center",
            )
            ax.axhline(0, color="#BDBDBD", linewidth=0.8, zorder=0)
            ax.set_xlabel("")
            ax.set_ylabel(
                f"{module_name}\nactivity coefficient" if col == 0 else "",
                color=module_color if col == 0 else "#222222",
                fontweight="bold" if col == 0 else "normal",
            )
            ax.spines[["top", "right"]].set_visible(False)
            ax.tick_params(axis="x", length=0)
            if row == 0:
                ax.set_title(PROJECT_LABELS[project], pad=8, fontweight="bold")
    fig.suptitle(
        "Lineage-associated iModulons retain parental-cell-line differences within individual projects",
        fontsize=11.5,
        fontweight="bold",
        y=0.995,
    )
    fig.subplots_adjust(left=0.09, right=0.995, top=0.89, bottom=0.11)
    png = OUT / "Supplementary_Figure_2_multilineage_projects.png"
    pdf = OUT / "Supplementary_Figure_2_multilineage_projects.pdf"
    fig.savefig(png, dpi=400, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    summary = (
        data.groupby(["Project", "Cell line"], observed=True)
        .agg(n=(52, "size"), ECM1_mean=(52, "mean"), ECM2_mean=(65, "mean"))
        .reset_index()
    )
    summary.to_csv(
        OUT / "Supplementary_Figure_2_multilineage_projects_summary.csv", index=False
    )
    print(png)
    print(pdf)
    print(summary.to_string(index=False))


if __name__ == "__main__":
    draw()
