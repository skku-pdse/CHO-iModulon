"""Reproduce discovery genomic plots from the original coordinate snapshot."""

from pathlib import Path
import argparse
import json
import pandas as pd
import matplotlib

matplotlib.use("Agg")
from genomic_style import add_genome_offsets, plot_one_module

from figure2_paths import DATA, RESULTS
ROOT = DATA.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--imodulons", nargs="+", type=int, help="Default: manuscript examples iM8 and iM25."
    )
    args = parser.parse_args()
    data = DATA
    out = RESULTS / "lineage"
    out.mkdir(parents=True, exist_ok=True)
    weights = pd.read_csv(ROOT / "discovery iModulon/data/discovery_gene_weight_matrix.csv")
    weights.columns = ["gene"] + [int(float(c)) for c in weights.columns[1:]]
    members = pd.read_csv(ROOT / "discovery iModulon/data/discovery_gene_membership.csv")
    positions = pd.read_csv(data / "genomic_positions.csv.gz")
    order = json.loads((data / "chromosome_order.json").read_text())
    merged = positions.merge(weights, on="gene", how="inner")
    merged, offsets = add_genome_offsets(merged, order)
    for module in args.imodulons or [8, 25]:
        frame = merged.copy()
        frame["highlight"] = frame.gene.isin(
            set(members.loc[members.iModulon.eq(module), "Gene"])
        )
        plot_one_module(
            frame,
            offsets,
            module,
            out,
            highlight_col="highlight",
            highlight_label=f"iModulon {module} genes",
        )
        print(f"Completed iM{module}", flush=True)


if __name__ == "__main__":
    main()
