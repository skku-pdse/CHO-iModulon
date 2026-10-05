"""Project external TPM or preprocessed expression onto the discovery gene-weight basis."""
from pathlib import Path
import argparse
import re
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler


def sample_key(project, sample):
    project = str(project).strip().replace("\u200b", "")
    project = {"R17704": "NISTCHO_Characterisation", "Project_BOKU": "NISTCHO_Characterisation"}.get(project, project)
    sample = re.sub(r"\.\d+$", "", str(sample).strip().replace("\u200b", ""))
    if project == "2020_02_12_ZeN" and sample == "D2-Z-ZeLa-2_S17":
        sample = "D3-Z-ZeLa-2_S17"
    return project, sample


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("expression", type=Path, nargs="?", help="Preprocessed genes-by-samples CSV")
    parser.add_argument("--tpm", type=Path, help="TPM TSV(.gz), with project/sample header rows")
    parser.add_argument("--samples", type=Path, help="Select samples using Project and SampleID columns in a metadata CSV")
    parser.add_argument("--basis", type=Path, default=Path(__file__).resolve().parents[1] / "discovery iModulon/data/discovery_gene_weight_matrix.csv")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "projected_activity.csv")
    args = parser.parse_args()
    if (args.expression is None) == (args.tpm is None):
        parser.error("Provide either a preprocessed expression CSV or --tpm.")
    if args.samples and not args.tpm:
        parser.error("--samples requires --tpm")
    basis = pd.read_csv(args.basis, index_col=0)
    if args.tpm:
        expression = pd.read_csv(args.tpm, sep="\t", index_col=0, header=[0, 1])
        if args.samples:
            expression.columns = pd.MultiIndex.from_tuples([sample_key(*c) for c in expression.columns])
            metadata = pd.read_csv(args.samples)
            keys = [sample_key(p, s) for p, s in zip(metadata.Project, metadata.SampleID)]
            if len(keys) != len(set(keys)):
                raise ValueError("Duplicate sample keys in metadata")
            expression = expression.loc[:, keys]
    else:
        expression = pd.read_csv(args.expression, index_col=0)
    assert basis.index.is_unique and expression.index.is_unique
    assert basis.columns.is_unique and expression.columns.is_unique
    missing = basis.index.difference(expression.index)
    if len(missing):
        raise ValueError(f"Expression is missing {len(missing)} basis genes")
    values = expression.loc[basis.index].to_numpy(float)
    assert np.isfinite(values).all() and np.isfinite(basis.to_numpy()).all()
    if args.tpm:
        if (values < 0).any():
            raise ValueError("TPM must be nonnegative")
        values = StandardScaler().fit_transform(np.log1p(values))
    activity, _, rank, _ = np.linalg.lstsq(basis.to_numpy(float), values, rcond=None)
    if rank < basis.shape[1]:
        raise ValueError("Basis is not full column rank")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(activity, index=basis.columns, columns=expression.columns).rename_axis("iModulon").to_csv(args.output)
    print(f"Saved {activity.shape[0]} components x {activity.shape[1]} samples")


if __name__ == "__main__":
    main()
