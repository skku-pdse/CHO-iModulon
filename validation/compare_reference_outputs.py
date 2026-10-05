"""Compare recalculated outputs with the supplied publication tables."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from projection.project import sample_key

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "outputs")
    args = parser.parse_args()
    rows = []

    def compare_table(actual, reference):
        if not actual.exists():
            return
        try:
            pd.testing.assert_frame_equal(pd.read_csv(actual), pd.read_csv(reference),
                                          check_dtype=False, check_exact=False, rtol=1e-9, atol=1e-10)
            status, detail = "match", ""
        except AssertionError as error:
            status, detail = "different", str(error).splitlines()[0]
        rows.append({"output": actual.name, "status": status, "detail": detail})

    for cohort in ("discovery", "expanded"):
        for folder, analysis in [("genomic", "chromosome_enrichment"), ("conditional association", "metadata_association")]:
            name = f"{cohort}_{analysis}_results.csv"
            compare_table(args.output_dir / cohort / name, ROOT / f"{cohort} iModulon" / folder / name)
        path = args.output_dir / cohort / f"{cohort}_membership_genes.csv"
        if path.exists():
            actual = pd.read_csv(path)
            expected = pd.read_csv(ROOT / f"{cohort} iModulon/data/{cohort}_gene_membership.csv")
            matches = set(zip(actual.iModulon, actual.Gene)) == set(zip(expected.iModulon, expected.Gene))
            rows.append({"output": path.name, "status": "match" if matches else "different", "detail": "gene-component pairs"})
    for path in (args.output_dir / "comparison").glob("*.csv"):
        reference = ROOT / "expanded iModulon/iModulon comparison" / path.name
        if reference.exists():
            compare_table(path, reference)
    path = args.output_dir / "projection/reproduced_projected_activities.csv"
    if path.exists():
        actual = pd.read_csv(path, index_col=0, header=[0, 1])
        expected = pd.read_csv(ROOT / "projection/data/projected_activity_matrix.csv", index_col=0, header=[0, 1])
        for frame in (actual, expected):
            frame.columns = pd.MultiIndex.from_tuples([sample_key(*c) for c in frame.columns])
        actual = actual.loc[expected.index, expected.columns]
        different = ~np.isclose(actual, expected, rtol=1e-8, atol=1e-8).all(axis=0)
        details = []
        for column in np.flatnonzero(different):
            details.append({"Project": expected.columns[column][0], "SampleID": expected.columns[column][1],
                            "max_absolute_difference": float(abs(actual.iloc[:, column]-expected.iloc[:, column]).max())})
        difference_path = path.parent / "projection_reference_differences.csv"
        if details:
            pd.DataFrame(details).to_csv(difference_path, index=False)
        elif difference_path.exists():
            difference_path.unlink()
        rows.append({"output": path.name, "status": "different" if different.any() else "match",
                     "detail": f"{different.sum()} of {len(different)} sample columns differ"})
    table = pd.DataFrame(rows, columns=["output", "status", "detail"])
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_dir / "reference_comparison.csv", index=False)
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()
