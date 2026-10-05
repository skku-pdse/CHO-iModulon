"""Highest version: revised validation workbook, 377 retained samples.

Adapted from figure4D_projected_temperature.py. Original inputs remain unchanged.
"""

import matplotlib

matplotlib.use("Agg")
from figure4_paths import DATA, RESULTS
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = DATA
XLSX = DATA / "validation.xlsx"
OUT = RESULTS / "external_projection"
OUT.mkdir(parents=True, exist_ok=True)
META, ACTIVITY, IM = ("Supplementary Table 14", "Supplementary Table 15", 39)
TS_PROJECTS = ("PRJNA778050", "PRJEB33024", "PRJNA1130622", "Project_BOKU")
EXCLUDE_PROJECTS = ("PRJNA214684",)
PROJECT_COLORS = {
    "PRJNA778050": "#c8553d",
    "PRJEB33024": "#e69f00",
    "PRJNA1130622": "#009E73",
    "Project_BOKU": "#7E57C2",
}
OTHER_COLOR = "#9aa3ad"
LABEL_FONT, TICK_FONT, LEGEND_FONT = (14, 10, 14)


def sample_id(value):
    s = str(value).strip()
    return s[:-2] if s.endswith(".0") and s[:-2].isdigit() else s


def load_data():
    meta = pd.read_excel(XLSX, sheet_name=META)
    meta.columns = meta.columns.astype(str).str.strip()
    meta["SampleID"] = meta["SampleID"].map(sample_id)
    raw = pd.read_excel(XLSX, sheet_name=ACTIVITY)
    ids = raw.iloc[0, 1:].map(sample_id)
    activity = raw.iloc[1:].copy()
    activity.iloc[:, 0] = pd.to_numeric(activity.iloc[:, 0], errors="coerce")
    row = activity.loc[activity.iloc[:, 0].eq(IM)]
    if len(row) != 1:
        raise ValueError(f"Expected exactly one iM{IM} row.")
    values = pd.to_numeric(row.iloc[0, 1:], errors="coerce").to_numpy()
    data = pd.DataFrame({"SampleID": ids, "Activity": values}).merge(
        meta, on="SampleID", how="left", validate="one_to_one"
    )
    if data.Project.isna().any():
        missing = data.loc[data.Project.isna(), "SampleID"].tolist()
        raise ValueError(
            f"{len(missing)} activity samples absent from {META}: {missing}"
        )
    dropped = data.Project.isin(EXCLUDE_PROJECTS).sum()
    data = data.loc[~data.Project.isin(EXCLUDE_PROJECTS)].copy()
    print(
        f"[info] loaded {len(data) + dropped} samples; excluded {dropped} ({', '.join(EXCLUDE_PROJECTS)}); {len(data)} retained."
    )
    return data


def plot(data):
    ts_project = data.Project.isin(TS_PROJECTS)
    shifted = data["Temperature shift"].eq("YES")
    groups = [
        ("Control\n(other projects)", data.loc[~ts_project]),
        ("Control\n(matched)", data.loc[ts_project & ~shifted]),
        ("Hypothermia\n", data.loc[ts_project & shifted]),
    ]
    fig, ax = plt.subplots(figsize=(5.0, 4.6))
    rng = np.random.default_rng(7)
    for pos, (label, subset) in enumerate(groups):
        subset = subset.dropna(subset=["Activity"])
        values = subset.Activity.to_numpy()
        q1, median, q3 = np.percentile(values, [25, 50, 75])
        print(
            f"[n] {label.strip().replace(chr(10), ' ')}: n={len(values)}  median={median:.3f}  IQR=[{q1:.3f}, {q3:.3f}]"
        )
        ax.add_patch(
            plt.Rectangle(
                (pos - 0.22, q1),
                0.44,
                q3 - q1,
                fill=True,
                facecolor="white",
                edgecolor="#3b3b3b",
                linewidth=0.7,
                alpha=0.9,
                zorder=2,
            )
        )
        ax.plot(
            [pos - 0.22, pos + 0.22],
            [median, median],
            color="#3b3b3b",
            linewidth=1.2,
            zorder=3,
        )
        x = pos + rng.uniform(
            -0.18 if pos == 0 else -0.12, 0.18 if pos == 0 else 0.12, len(values)
        )
        colors = subset.Project.map(PROJECT_COLORS).fillna(OTHER_COLOR)
        boku = subset.Project.eq("Project_BOKU")
        if pos == 0:
            ax.scatter(
                x,
                values,
                s=5,
                c=OTHER_COLOR,
                alpha=0.35,
                edgecolor="none",
                linewidth=0,
                zorder=4,
            )
        else:
            ax.scatter(
                x[boku],
                values[boku],
                s=42,
                c=colors[boku],
                alpha=0.92,
                edgecolor="white",
                linewidth=0.5,
                zorder=4,
            )
            ax.scatter(
                x[~boku],
                values[~boku],
                s=42,
                c=colors[~boku],
                alpha=0.92,
                edgecolor="white",
                linewidth=0.5,
                zorder=4,
            )
        groups[pos] = (f"{label}\n(n={len(values)})", subset)
    focus = data.loc[ts_project, "Activity"].dropna()
    pad = 0.1 * (focus.max() - focus.min())
    ax.set_ylim(focus.min() - pad, focus.max() + 0.6 * pad)
    ax.axhline(0, color="#666", linewidth=0.7, linestyle="--", alpha=0.6, zorder=1)
    ax.set_xticks(range(3))
    ax.set_xticklabels([label for label, _ in groups], fontsize=TICK_FONT)
    ax.set_ylabel("iM39 activity (A')", fontsize=LABEL_FONT)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="y", labelsize=TICK_FONT)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=PROJECT_COLORS[p],
            markeredgecolor="none",
            markersize=12,
            label=p,
        )
        for p in ("PRJNA778050", "PRJEB33024", "PRJNA1130622", "Project_BOKU")
    ]
    handles.append(
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=OTHER_COLOR,
            markeredgecolor="none",
            markersize=12,
            label="Other projects",
        )
    )
    ax.legend(
        handles=handles,
        loc="upper left",
        frameon=False,
        fontsize=LEGEND_FONT,
        handletextpad=0.45,
        labelspacing=0.28,
        borderaxespad=0.55,
    )
    fig.tight_layout()
    ax.set_title("Temperature Downshift (iM39)")
    fig.savefig(OUT / "figure4D_projected_temperature.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] saved -> {OUT / 'figure4D_projected_temperature.png'}")


if __name__ == "__main__":
    plot(load_data())
