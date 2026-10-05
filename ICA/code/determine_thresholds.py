"""Gene membership from the last large gap in sorted absolute ICA weights."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def determine_thresholds(matrix, gap_ratio=0.01):
    """Return long-format members and per-component threshold diagnostics.

    Reproduces the initial gap-method cells of
    imodulon_threshold_determination.ipynb, including >= comparisons,
    original gene order, and minimum ranks for ties.
    """
    if matrix.empty or not matrix.index.is_unique or not matrix.columns.is_unique:
        raise ValueError("Matrix must be nonempty with unique genes and components.")
    if not np.isfinite(matrix.to_numpy(dtype=float)).all():
        raise ValueError("All matrix coefficients must be finite numbers.")
    if not np.isfinite(gap_ratio) or gap_ratio <= 0:
        raise ValueError("gap_ratio must be positive.")
    members, diagnostics = [], []
    for component in matrix.columns:
        weights = matrix[component]
        absolute = weights.abs()
        values = absolute.sort_values(ascending=False).to_numpy()
        gaps = -np.diff(values)
        gap_threshold = values[0] * gap_ratio
        hits = np.flatnonzero(gaps >= gap_threshold)
        cutoff = int(hits[-1]) if len(hits) else None
        threshold = values[cutoff] if cutoff is not None else np.nan
        selected = weights[absolute >= threshold] if cutoff is not None else weights.iloc[:0]
        members.append(pd.DataFrame({
            "iModulon": component,
            "Gene": selected.index,
            "Rank": selected.abs().rank(ascending=False, method="min").to_numpy(),
            "Loading": selected.to_numpy(),
        }))
        diagnostics.append({
            "iModulon": component, "max_abs_coefficient": values[0],
            "gap_ratio": gap_ratio,
            "gap_threshold": gap_threshold,
            "cutoff_rank": cutoff + 1 if cutoff is not None else np.nan,
            "abs_coefficient_threshold": threshold,
            "selected_gap": gaps[cutoff] if cutoff is not None else np.nan,
            "n_genes": len(selected),
        })
    return pd.concat(members, ignore_index=True), pd.DataFrame(diagnostics)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("matrix", type=Path, help="M.csv: genes by components")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--prefix", default="membership", help="Unique output filename prefix")
    parser.add_argument("--gap-ratio", type=float, default=0.01)
    args = parser.parse_args()
    matrix = pd.read_csv(args.matrix, index_col=0)
    members, thresholds = determine_thresholds(matrix, args.gap_ratio)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    members.to_csv(args.output_dir / f"{args.prefix}_genes.csv", index=False)
    thresholds.to_csv(args.output_dir / f"{args.prefix}_thresholds.csv", index=False)
    print(f"{len(thresholds)} components; {len(members)} membership rows")


if __name__ == "__main__":
    main()
