"""Recompute comparisons; retain five publication tables and no intermediate files."""
from pathlib import Path
import argparse
import os
import subprocess
import sys
import tempfile
import shutil
import pandas as pd

HERE = Path(__file__).resolve().parent


def calculate(work):
    env = dict(os.environ, IMODULON_COMPARISON_WORK=str(work), PYTHONIOENCODING="utf-8")
    for script in ["compare_components.py", "analyze_activities.py"]:
        subprocess.run([sys.executable, str(HERE / script)], env=env, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=HERE)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="imodulon_comparison_") as directory:
        work = Path(directory)
        calculate(work)
        membership = pd.read_csv(work / "all_pairwise_membership.csv")
        weights = pd.read_csv(work / "pairwise_gene_weight_correlations.csv", index_col=0)
        weights.index = weights.index.astype(int)
        weights.columns = weights.columns.astype(int)
        membership["weight_r"] = [weights.loc[d, e] for d, e in zip(membership.discovery_iModulon, membership.expanded_iModulon)]
        membership.to_csv(args.output_dir / "component_comparison.csv", index=False)
        for source, target in {
            "activity_and_membership_comparison.csv": "activity_comparison.csv",
            "external_projection_correlations.csv": "external_activity.csv",
            "DXB11_exact_permutation.csv": "DXB11_permutation.csv",
            "representative_member_gene_changes.csv": "member_changes.csv",
        }.items():
            shutil.copy2(work / source, args.output_dir / target)
    print("Saved five tables; intermediate files removed automatically.")


if __name__ == "__main__":
    main()
