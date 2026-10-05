"""Current revised activity heatmap.

Adapted from figure3F_representative_activity_heatmap.py. Original inputs remain unchanged.
"""

import matplotlib

matplotlib.use("Agg")
from figure1_paths import DATA, RESULTS
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Rectangle

ROOT = DATA
IMOD_XLSX = DATA / "discovery.xlsx"
META_XLSX = DATA / "metadata.xlsx"
OUT_DIR = RESULTS / "overview"
OUT_DIR.mkdir(parents=True, exist_ok=True)
ROWS = [
    ("Functional / regulatory", 4, "UPR iM4"),
    ("Functional / regulatory", 9, "Translation-2 iM9"),
    ("Functional / regulatory", 16, "Cell Cycle-1 iM16"),
    ("Functional / regulatory", 2, "Lysosome-1 iM2"),
    ("Lineage-associated", 52, "ECM-1 iM52"),
    ("Lineage-associated", 65, "ECM-2 iM65"),
    ("Genomic locus", 8, "Ephrin Locus iM8"),
    ("Genomic locus", 25, "Histone iM25"),
    ("Process-responsive", 39, "TDS iM39"),
    ("Process-responsive", 87, "ISG-1 iM87"),
    ("Process-responsive", 54, "ATF4 iM54"),
]
ROW_BLOCK_ORDER = [
    "Functional / regulatory",
    "Lineage-associated",
    "Genomic locus",
    "Process-responsive",
]


def fda_class(cond):
    s = str(cond)
    has_V = any(
        (tok in s for tok in ["Vm", "Vmedia", "Vifn", "VIFN", "Vpoly", "VpolyIC"])
    )
    has_im = "ifn" in s.lower() or "poly" in s.lower()
    if has_V and has_im:
        return "Virus+immune"
    if has_V and (not has_im):
        return "Virus only"
    if not has_V and has_im:
        return "IFN/polyIC only"
    return "FDA control"


def in_fda(md):
    return md["Project"].astype(str).str.contains("FDA", case=False, na=False)


COLUMNS = [
    ("Lineage", "CHO-K1", lambda md: md["Cell line"] == "CHO-K1"),
    ("Lineage", "CHO-S", lambda md: md["Cell line"] == "CHO-S"),
    ("Lineage", "CHO-DG44", lambda md: md["Cell line"] == "CHO-DG44"),
    ("Lineage", "CHO-ZeLa", lambda md: md["Cell line"] == "CHO-ZeLa"),
    ("Culture phase", "Exponential", lambda md: md["Phase"] == "E"),
    ("Culture phase", "Stationary", lambda md: md["Phase"].isin(["S"])),
    ("Culture phase", "Decline", lambda md: md["Phase"].isin(["D"])),
    ("Temperature", "NTS", lambda md: md["Temperature shift"] == "NO"),
    ("Temperature", "TS", lambda md: md["Temperature shift"] == "YES"),
    (
        "Antiviral",
        "FDA ctrl",
        lambda md: in_fda(md) & (md["Condition"].apply(fda_class) == "FDA control"),
    ),
    (
        "Antiviral",
        "IFN/polyIC",
        lambda md: in_fda(md) & (md["Condition"].apply(fda_class) == "IFN/polyIC only"),
    ),
    (
        "Antiviral",
        "Virus",
        lambda md: in_fda(md) & (md["Condition"].apply(fda_class) == "Virus only"),
    ),
    (
        "Antiviral",
        "Virus+immune",
        lambda md: in_fda(md) & (md["Condition"].apply(fda_class) == "Virus+immune"),
    ),
]
COL_BLOCK_ORDER = ["Lineage", "Culture phase", "Temperature", "Antiviral"]
COL_BLOCK_COLOR = {
    "Lineage": "#0072B2",
    "Culture phase": "#CC79A7",
    "Temperature": "#E69F00",
    "Antiviral": "#D55E00",
}
print("[1/4] Loading A and metadata ...")
A_re = pd.read_excel(IMOD_XLSX, sheet_name="A_reorganized")
SID_COL = "Unnamed: 1"
imod_cols = [c for c in A_re.columns if isinstance(c, (int, np.integer))]
A_re = A_re.dropna(subset=[SID_COL]).drop_duplicates(subset=SID_COL, keep="first")
A_re[SID_COL] = A_re[SID_COL].astype(str)
A = A_re.set_index(SID_COL)[imod_cols].apply(pd.to_numeric, errors="coerce").T
A.index = A.index.astype(int)
md = pd.read_excel(META_XLSX, sheet_name="Main_v3")
md.columns = [c.strip() for c in md.columns]
md["SampleID"] = md["SampleID"].astype(str)
md = md.drop_duplicates(subset="SampleID", keep="first").set_index("SampleID")
common = [s for s in A.columns if s in md.index]
A = A[common]
md = md.loc[common]
print(f"     n_samples aligned = {len(common)}")
A_z = A.sub(A.mean(axis=1), axis=0).div(A.std(axis=1, ddof=1), axis=0)
print("[2/4] Computing group means ...")
n_per_col = []
group_means = []
col_labels = []
col_blocks = []
for block, label, sel in COLUMNS:
    mask = sel(md).fillna(False).values
    samples = md.index[mask]
    n = len(samples)
    n_per_col.append(n)
    col_labels.append(label)
    col_blocks.append(block)
    if n == 0:
        group_means.append(pd.Series(np.nan, index=A_z.index))
    else:
        group_means.append(A_z[samples].mean(axis=1))
    print(f"     {block:>14s} | {label:<14s}  n={n:>3d}")
H_full = pd.concat(group_means, axis=1)
H_full.columns = col_labels
row_imods = [r[1] for r in ROWS]
row_labels = [r[2] for r in ROWS]
H = H_full.loc[row_imods]
H.index = row_labels
assert A.columns.equals(md.index), "Activity and metadata order mismatch"
assert len(set(row_imods)) == len(row_imods)
assert np.isfinite(H.to_numpy()).all(), "Missing or non-finite group means"
assert all((n > 0 for n in n_per_col)), "Empty sample group"
out_csv = OUT_DIR / "figure1G_heatmap_groupmean_zscore.csv"
H.to_csv(out_csv)
print(f"     wrote {out_csv}")
print("[3/4] Plotting ...")
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 35,
        "axes.labelsize": 40,
        "xtick.labelsize": 32.5,
        "ytick.labelsize": 32.5,
        "legend.fontsize": 35,
        "axes.linewidth": 3.0,
        "xtick.major.width": 3.0,
        "ytick.major.width": 3.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)
n_rows, n_cols = H.shape
fig, ax = plt.subplots(figsize=(35.0, 20.0), dpi=150)
vmax = float(np.nanpercentile(np.abs(H.values), 98))
vmax = max(vmax, 1.0)
vmin = -vmax
cmap = plt.get_cmap("RdBu_r").copy()
cmap.set_bad(color="#dddddd")
im = ax.imshow(H.values, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
for i in range(n_rows):
    for j in range(n_cols):
        v = H.values[i, j]
        if np.isnan(v):
            ax.text(
                j, i, "n/a", ha="center", va="center", fontsize=27.5, color="#666666"
            )
        elif abs(v) >= 1:
            ax.text(
                j,
                i,
                f"{v:+.1f}",
                ha="center",
                va="center",
                fontsize=30,
                color="white" if abs(v) > 1.6 else "#222222",
            )
ax.set_xticks(np.arange(n_cols))
ax.set_xticklabels(
    [f"{lab}\n(n={n})" for lab, n in zip(col_labels, n_per_col)],
    rotation=45,
    ha="right",
    fontsize=32.5,
)
ax.set_yticks(np.arange(n_rows))
ax.set_yticklabels(H.index, fontsize=40)
ax.tick_params(axis="both", length=10, pad=7.5)
prev_block = col_blocks[0]
for j in range(1, n_cols):
    if col_blocks[j] != prev_block:
        ax.axvline(j - 0.5, color="white", lw=10.0)
        ax.axvline(j - 0.5, color="#222222", lw=2.5)
    prev_block = col_blocks[j]
row_blocks = [r[0] for r in ROWS]
prev_block = row_blocks[0]
for i in range(1, n_rows):
    if row_blocks[i] != prev_block:
        ax.axhline(i - 0.5, color="white", lw=10.0)
        ax.axhline(i - 0.5, color="#222222", lw=2.5)
    prev_block = row_blocks[i]


def block_label_positions(blocks):
    out = []
    start = 0
    cur = blocks[0]
    for j in range(1, len(blocks) + 1):
        if j == len(blocks) or blocks[j] != cur:
            out.append((cur, start, j - 1))
            if j < len(blocks):
                cur = blocks[j]
                start = j
    return out


for block_name, j0, j1 in block_label_positions(col_blocks):
    cx = (j0 + j1) / 2
    color = COL_BLOCK_COLOR.get(block_name, "#444444")
    ax.text(
        cx,
        -1.4,
        block_name,
        ha="center",
        va="bottom",
        fontsize=40,
        fontweight="bold",
        color=color,
        transform=ax.transData,
    )
    ax.plot([j0 - 0.45, j1 + 0.45], [-0.85, -0.85], color=color, lw=10.0, clip_on=False)
cbar = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.012, extend="both")
cbar.set_label("Group-mean row-z (activity)", fontsize=35)
cbar.ax.tick_params(labelsize=32.5, length=10)
plt.subplots_adjust(left=0.22, right=0.97, top=0.84, bottom=0.22)
png = OUT_DIR / "figure1G_representative_activity_heatmap.png"
pdf = OUT_DIR / "figure1G_representative_activity_heatmap.pdf"
fig.savefig(png, dpi=150, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
plt.close(fig)
print(f"     wrote {png}")
print("\n[4/4] Sanity check on expected patterns:")


def show(imod_id, cols):
    if imod_id not in row_imods:
        return
    lab = ROWS[row_imods.index(imod_id)][2]
    vals = H_full.loc[imod_id, cols]
    print(f"  {lab}")
    for c, v in vals.items():
        print(f"     {c:<14s} z={v:+.2f}")


show(39, ["NTS", "TS"])
show(87, ["FDA ctrl", "IFN/polyIC", "Virus", "Virus+immune"])
show(52, ["CHO-K1", "CHO-S", "CHO-DG44", "CHO-ZeLa"])
show(65, ["CHO-K1", "CHO-S", "CHO-DG44", "CHO-ZeLa"])
show(4, ["Exponential", "Stationary", "Decline", "NTS", "TS"])
show(57, ["Exponential", "Stationary", "Decline", "CHO-DG44"])
show(5, col_labels)
show(25, col_labels)
print("\nDone.")
