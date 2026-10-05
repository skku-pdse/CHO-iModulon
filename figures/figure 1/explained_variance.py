"""Plot archived single-component reconstruction EV; proxy variant excluded.

Adapted from figure1D_explained_variance.py. Original inputs remain unchanged.
"""

import matplotlib

matplotlib.use("Agg")
from figure1_paths import DATA, RESULTS
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# Plot the manuscript's archived reconstruction-EV values. ICA preprocessing
# and EV recalculation are separate from this figure-only script.
OUT_DIR = RESULTS / "overview"
OUT_DIR.mkdir(parents=True, exist_ok=True)
ev_df = pd.read_csv(DATA / "explained_variance_reference.csv")
ev_df.columns = ev_df.columns.str.strip()
ev_df = ev_df.rename(columns={"imodulon": "iModulon"})
ev_df = ev_df.sort_values("explained_variance", ascending=False).reset_index(drop=True)
ev_df["cumulative"] = ev_df.explained_variance.cumsum()

per_var = ev_df["explained_variance"].values.astype(float)
cum_var = per_var.cumsum()
n = len(per_var)
x = np.arange(1, n + 1)
fig, ax1 = plt.subplots(figsize=(7.2, 6.0), dpi=300)
ax1.bar(x, per_var, width=0.82, color="#c7d8ee", edgecolor="none", zorder=1)
ax1.set_yscale("log")
positive = per_var[per_var > 0]
ymin = max(positive.min() / 2, 1e-06) if len(positive) else 1e-06
ymax = max(per_var.max() * 2, 0.001)
ax1.set_ylim(ymin, ymax)
ax1.set_xlim(0.25, n + 0.75)
ax1.set_xlabel("iModulons ranked by EV")
ax1.set_ylabel("Per-iModulon EV (%, log)", color="#444444")
ax1.tick_params(axis="y", colors="#444444", pad=7.5)
ax1.tick_params(axis="x", pad=7.5)
ax1.spines["top"].set_visible(False)
ax1.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
CUM_COLOR = "#1f3b73"
ax2 = ax1.twinx()
ax2.plot(x, cum_var, color=CUM_COLOR, linewidth=5.0, zorder=3)
ax2.set_ylabel("Cumulative EV (%)", color=CUM_COLOR)
ax2.tick_params(axis="y", colors=CUM_COLOR, pad=7.5, length=12.5, width=3.0)
ax2.set_ylim(0, max(100, cum_var.max() * 1.05))
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_color(CUM_COLOR)
ax2.spines["right"].set_linewidth(3.0)
plt.tight_layout(pad=2.0)
png = OUT_DIR / "figure1D_explained_variance.png"
pdf = OUT_DIR / "figure1D_explained_variance.pdf"
fig.savefig(png, dpi=300, bbox_inches="tight")
fig.savefig(pdf, bbox_inches="tight")
plt.close(fig)
print(f"Wrote {png}")
print(f"Wrote {pdf}")
