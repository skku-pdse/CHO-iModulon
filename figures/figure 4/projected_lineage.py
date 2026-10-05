"""Highest version: revised validation workbook and NIST-CHO.

Adapted from figure4BC_projected_lineage.py. Original inputs remain unchanged.
"""

import matplotlib

matplotlib.use("Agg")
from figure4_paths import DATA, RESULTS
from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

ROOT = DATA
XLSX = DATA / "validation.xlsx"
OUT = RESULTS / "external_projection"
OUT.mkdir(parents=True, exist_ok=True)
META, ACTIVITY = ("Supplementary Table 14", "Supplementary Table 15")
EXCLUDE_PROJECTS = ("PRJNA214684",)
LINEAGES = ["K1", "DG44", "S", "DXB11", "NIST-CHO"]
LABELS = {
    "K1": "CHO-K1",
    "DG44": "CHO-DG44",
    "S": "CHO-S",
    "DXB11": "CHO-DXB11",
    "NIST-CHO": "NIST-CHO",
}
PALETTE = {
    52: {"violin": "#E89A8D", "dot": "#c0392b"},
    65: {"violin": "#8FB6E0", "dot": "#2f69b0"},
}
LABEL_FONT, TICK_FONT = (14, 10)


def sample_id(value):
    s = str(value).strip()
    return s[:-2] if s.endswith(".0") and s[:-2].isdigit() else s


def lineage(value):
    s = str(value).strip()
    return {
        "CHO-K1": "K1",
        "K1": "K1",
        "CHO-DG44": "DG44",
        "DG44": "DG44",
        "CHO-S": "S",
        "S": "S",
        "CHO-DXB11": "DXB11",
        "DXB11": "DXB11",
        "NIST-CHO": "NIST-CHO",
    }.get(s)


def three_line_label(key, n):
    """Use a shared three-line x-axis format to prevent lineage-label overlap."""
    prefix, suffix = LABELS[key].split("-", 1)
    if suffix == "S" or suffix == "K1":
        return f"{prefix}-{suffix}\n\n(n={n})"
    else:
        return f"{prefix}-\n{suffix}\n(n={n})"


def load_data(im):
    meta = pd.read_excel(XLSX, sheet_name=META)
    meta.columns = meta.columns.astype(str).str.strip()
    meta.SampleID = meta.SampleID.map(sample_id)
    raw = pd.read_excel(XLSX, sheet_name=ACTIVITY)
    ids = raw.iloc[0, 1:].map(sample_id)
    activity = raw.iloc[1:].copy()
    activity.iloc[:, 0] = pd.to_numeric(activity.iloc[:, 0], errors="coerce")
    row = activity.loc[activity.iloc[:, 0].eq(im)]
    if len(row) != 1:
        raise ValueError(f"Expected exactly one iM{im} row.")
    df = pd.DataFrame(
        {
            "SampleID": ids,
            "Activity": pd.to_numeric(row.iloc[0, 1:], errors="coerce").to_numpy(),
        }
    ).merge(
        meta[["SampleID", "Project", "Cell line"]],
        on="SampleID",
        how="left",
        validate="one_to_one",
    )
    if df.Project.isna().any():
        missing = df.loc[df.Project.isna(), "SampleID"].tolist()
        raise ValueError(
            f"{len(missing)} activity samples absent from {META}: {missing}"
        )
    dropped = df.Project.isin(EXCLUDE_PROJECTS).sum()
    df = df.loc[~df.Project.isin(EXCLUDE_PROJECTS)].copy()
    print(
        f"[info] iM{im}: {len(df) + dropped} samples loaded, {dropped} excluded ({', '.join(EXCLUDE_PROJECTS)}), {len(df)} retained."
    )
    df["lineage"] = df["Cell line"].map(lineage)
    return df.dropna(subset=["Activity", "lineage"])


def draw(ax, df, im, title):
    pal = PALETTE[im]
    for key in LINEAGES:
        sub = df.loc[df.lineage == key, "Activity"]
        print(f"[n] iM{im} {LABELS[key]}: n={len(sub)}  median={sub.median():.3f}")
    sns.violinplot(
        data=df,
        x="lineage",
        y="Activity",
        order=LINEAGES,
        inner=None,
        cut=0,
        linewidth=1.4,
        color=pal["violin"],
        saturation=1,
        ax=ax,
    )
    sns.stripplot(
        data=df,
        x="lineage",
        y="Activity",
        order=LINEAGES,
        color=pal["dot"],
        alpha=0.8,
        size=2.7,
        jitter=0.18,
        linewidth=0.5,
        ax=ax,
    )
    ax.axhline(0, color=".75", linewidth=1, zorder=0)
    ax.set_xticks(range(len(LINEAGES)))
    ax.set_xticklabels(
        [three_line_label(x, (df.lineage == x).sum()) for x in LINEAGES],
        fontsize=TICK_FONT,
    )
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("Activity (A')", fontsize=LABEL_FONT)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="y", labelsize=TICK_FONT)


if __name__ == "__main__":
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 6), gridspec_kw={"wspace": 0.28})
    draw(axes[0], load_data(52), 52, "ECM-1 (iM52)")
    draw(axes[1], load_data(65), 65, "ECM-2 (iM65)")
    fig.subplots_adjust(bottom=0.2, top=0.9, wspace=0.28)
    fig.savefig(
        OUT / "figure4BC_projected_lineage.png", dpi=300, bbox_inches="tight"
    )
    plt.close(fig)
    print(f"[OK] saved -> {OUT / 'figure4BC_projected_lineage.png'}")
