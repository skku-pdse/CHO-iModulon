"""Check the shapes, identifiers, and numerical validity of packaged analysis inputs."""

from pathlib import Path
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from projection.project import sample_key

ROOT = Path(__file__).resolve().parents[1]


def main():
    for cohort, modules, samples in [("discovery", 105, 778), ("expanded", 148, 1030)]:
        base = ROOT / f"{cohort} iModulon/data"
        m = pd.read_csv(base / f"{cohort}_gene_weight_matrix.csv", index_col=0)
        a = pd.read_csv(base / f"{cohort}_activity_matrix.csv", index_col=0, header=[0, 1])
        metadata = pd.read_csv(base / f"{cohort}_sample_metadata.csv")
        members = pd.read_csv(base / f"{cohort}_gene_membership.csv")
        if m.shape != (19812, modules) or a.shape != (modules, samples):
            raise ValueError(f"Unexpected {cohort} matrix dimensions")
        if not m.index.is_unique or not np.isfinite(m.to_numpy()).all() or not np.isfinite(a.to_numpy()).all():
            raise ValueError(f"Invalid {cohort} matrix")
        activity_keys = [sample_key(*key) for key in a.columns]
        metadata_keys = [sample_key(p, s) for p, s in zip(metadata.Project, metadata.SampleID)]
        if len(set(activity_keys)) != samples or len(set(metadata_keys)) != samples or set(activity_keys) != set(metadata_keys):
            raise ValueError(f"{cohort} activity/metadata keys differ")
        if {int(float(v)) for v in m.columns} != {int(float(v)) for v in a.index}:
            raise ValueError(f"{cohort} component labels differ")
        if not set(members.Gene) <= set(m.index) or members.duplicated(["iModulon", "Gene"]).any():
            raise ValueError(f"Invalid {cohort} membership")
        print(f"{cohort}: {m.shape[0]} genes, {samples} samples, {modules} iModulons; inputs valid")
    projected = pd.read_csv(ROOT / "projection/data/projected_activity_matrix.csv", index_col=0, header=[0, 1])
    metadata = pd.read_csv(ROOT / "projection/data/projected_sample_metadata.csv")
    if projected.shape != (105, 255) or not np.isfinite(projected.to_numpy()).all():
        raise ValueError("Invalid projected activity matrix")
    keys = [sample_key(p, s) for p, s in zip(metadata.Project, metadata.SampleID)]
    if len(set(keys)) != 255 or set(keys) != {sample_key(*c) for c in projected.columns}:
        raise ValueError("Projected metadata keys differ")
    print("projection: 255 samples; inputs valid")


if __name__ == "__main__":
    main()
