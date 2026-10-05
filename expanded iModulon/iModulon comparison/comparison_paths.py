"""Resolve legacy reader interfaces to the canonical package data."""
from pathlib import Path
import os
import pandas as pd
ROOT = next(p for p in Path(__file__).resolve().parents if (p / "discovery iModulon").is_dir())
RESULTS = Path(os.environ["IMODULON_COMPARISON_WORK"])
RESULTS.mkdir(parents=True, exist_ok=True)
class Inputs:
    def __truediv__(self, name):
        if name == "metadata.csv": return ROOT / "expanded iModulon/data/expanded_sample_metadata.csv"
        if name == "projected_metadata.csv": return ROOT / "projection/data/projected_sample_metadata.csv"
        if name == "projected_A.csv.gz": return ROOT / "projection/data/projected_activity_matrix.csv"
        for cohort in ["discovery", "expanded"]:
            if name.startswith(cohort + "_"):
                tail = name[len(cohort)+1:].replace(".csv.gz", ".csv")
                return ROOT / f"{cohort} iModulon/data" / (cohort + "_" + {"A.csv": "activity_matrix.csv", "M.csv": "gene_weight_matrix.csv", "members.csv": "gene_membership.csv", "names.csv": "imodulon_names.csv"}[tail])
        raise KeyError(name)
DATA = Inputs()
def read_weights(cohort):
    matrix = pd.read_csv(ROOT / f"{cohort} iModulon/data/{cohort}_gene_weight_matrix.csv", index_col=0)
    matrix.columns = [int(float(c)) for c in matrix.columns]
    return matrix
