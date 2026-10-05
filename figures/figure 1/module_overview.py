"""Plot module sizes, category intersections, and EV-weighted functional families."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import squarify
from figure1_paths import DATA, RESULTS

table = pd.read_csv(DATA / "discovery_module_annotations.csv")
output = RESULTS / "module_overview"
output.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42})


def save(fig, name):
    for extension in ("png", "pdf"):
        fig.savefig(output / f"{name}.{extension}", dpi=300, bbox_inches="tight")
    plt.close(fig)


fig, ax = plt.subplots(figsize=(4.5, 3.2))
ax.hist(table["Total gene number"], bins=np.arange(0, 171, 10), color="#FF7B80", edgecolor="white")
ax.set(xlabel="Genes per iModulon", ylabel="Number of iModulons")
ax.spines[["top", "right"]].set_visible(False)
save(fig, "figure1B_module_sizes")

categories = ["Functional", "Conditional", "Regulatory", "RNA-biotype", "Uncharacterized", "Genomic"]
memberships = table.Category.map(lambda x: tuple(c.strip() for c in x.split(",")))
combinations = memberships.map(lambda x: tuple(c in x for c in categories)).value_counts()
combinations = sorted(combinations.items(), key=lambda pair: (sum(pair[0]), -pair[1], pair[0]))
fig = plt.figure(figsize=(8, 4.8))
grid = fig.add_gridspec(2, 2, width_ratios=[1.4, 5], height_ratios=[2, 1.5], hspace=.06, wspace=.05)
bars = fig.add_subplot(grid[0, 1]); dots = fig.add_subplot(grid[1, 1], sharex=bars)
totals = fig.add_subplot(grid[1, 0], sharey=dots)
x = np.arange(len(combinations))
counts = [count for _, count in combinations]
bars.bar(x, counts, color="black", width=.65)
for i, count in enumerate(counts): bars.text(i, count + .4, str(count), ha="center", fontsize=8)
bars.set_ylabel("Intersection size"); bars.tick_params(axis="x", bottom=False, labelbottom=False)
bars.set_ylim(0, max(counts) * 1.17)
for i, (combination, _) in enumerate(combinations):
    dots.scatter([i] * len(categories), np.arange(len(categories)), c="#D0D0D0", s=38)
    active = np.flatnonzero(combination)
    dots.plot([i] * len(active), active, "o-", color="black", markersize=6)
dots.set(yticks=np.arange(len(categories)), yticklabels=[], xticks=[], ylim=(-.5, len(categories)-.5))
total_counts = [sum(c in row for row in memberships) for c in categories]
totals.barh(np.arange(len(categories)), total_counts, color="black", height=.6)
totals.set_yticks(np.arange(len(categories)), categories); totals.invert_xaxis()
totals.set_xlabel("iModulons")
for axis in (bars, dots, totals): axis.spines[["top", "right"]].set_visible(False)
dots.spines[["left", "bottom"]].set_visible(False)
save(fig, "figure1C_category_intersections")

palette = {"Uncharacterized": "#D9D9D9", "ECM / extracellular": "#93C47D",
           "Stress / immune": "#6AA84F", "Genomic localization": "#76A5AF",
           "RNA-biotype": "#C00000", "Translation": "#9FC5E8", "Proteostasis": "#F8CBAD",
           "Cell cycle": "#FFD966", "Metabolism / redox": "#8E7CC3",
           "Regulatory-associated": "#C27BA0", "Chromatin / development": "#45818E",
           "Conditional / context": "#F6B26B", "Other functional": "#B4A7D6"}
families = table.groupby("broad_family").EV.sum().sort_values(ascending=False)
blocks = squarify.squarify(squarify.normalize_sizes(families.values, 100, 140), 0, 0, 100, 140)
fig, ax = plt.subplots(figsize=(6, 8))
for (family, mass), block in zip(families.items(), blocks):
    genes = table.loc[table.broad_family.eq(family)].sort_values("EV", ascending=False)
    tiles = squarify.squarify(squarify.normalize_sizes(genes.EV, block["dx"], block["dy"]),
                             block["x"], block["y"], block["dx"], block["dy"])
    for row, tile in zip(genes.itertuples(), tiles):
        ax.add_patch(Rectangle((tile["x"], tile["y"]), tile["dx"], tile["dy"],
                               facecolor=palette[family], edgecolor="white", linewidth=.5))
        if tile["dx"] >= 6 and tile["dy"] >= 4:
            ax.text(tile["x"]+tile["dx"]/2, tile["y"]+tile["dy"]/2, row.iModulon,
                    ha="center", va="center", fontsize=7)
    ax.add_patch(Rectangle((block["x"], block["y"]), block["dx"], block["dy"],
                           fill=False, edgecolor="white", linewidth=2))
ax.set(xlim=(0, 100), ylim=(0, 140)); ax.axis("off")
handles = [Rectangle((0, 0), 1, 1, color=palette[f]) for f in families.index]
fig.legend(handles, families.index, loc="lower center", ncol=2, fontsize=8, frameon=False)
fig.subplots_adjust(bottom=.18)
save(fig, "figure1E_functional_treemap")
