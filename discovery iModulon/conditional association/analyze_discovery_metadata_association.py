"""Nine metadata factors, omega squared, 1,000 permutations, global BH correction.

Writes one statistical table. Does not assign a Conditional category.
"""
from pathlib import Path
import argparse
import re
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
FIELDS = ["Temperature shift", "Time point (day)", "Engineered", "Phase", "Culture",
          "Producer", "Media", "Cell line", "Project"]


def sample_key(project, sample):
    project = {"R17704": "NISTCHO_Characterisation", "Project_BOKU": "NISTCHO_Characterisation"}.get(str(project).strip(), str(project).strip())
    sample = re.sub(r"\.\d+$", "", str(sample).strip())
    if project == "2020_02_12_ZeN" and sample == "D2-Z-ZeLa-2_S17":
        sample = "D3-Z-ZeLa-2_S17"
    return project, sample


def effect(values, labels):
    groups, index = np.unique(labels, return_inverse=True)
    n, k = values.shape[1], len(groups)
    if n < 6 or k < 2 or n <= k:
        return None
    counts = np.bincount(index).astype(float)
    means = np.column_stack([values[:, index == j].sum(axis=1) for j in range(k)]) / counts
    grand = values.mean(axis=1, keepdims=True)
    total = ((values - grand) ** 2).sum(axis=1)
    between = (counts * (means - grand) ** 2).sum(axis=1)
    within = np.maximum(total - between, 1e-30) / (n - k)
    f_stat = between / (k - 1) / within
    omega = np.maximum((between - (k - 1) * within) / (total + within), 0)
    return omega, f_stat, stats.f.sf(f_stat, k - 1, n - k), n, k


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "discovery_metadata_association_results.csv")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    raw = pd.read_csv(HERE.parent / "data/discovery_activity_matrix.csv", header=None)
    modules = pd.to_numeric(raw.iloc[2:, 0]).astype(int).to_numpy()
    activity = raw.iloc[2:, 1:].to_numpy(float)
    keys = [sample_key(p, s) for p, s in zip(raw.iloc[0, 1:], raw.iloc[1, 1:])]
    metadata = pd.read_csv(HERE.parent / "data/discovery_sample_metadata.csv")
    metadata.columns = metadata.columns.str.strip()
    for column in metadata.select_dtypes("object"):
        metadata[column] = metadata[column].map(lambda value: value.strip() if isinstance(value, str) else value).replace("", np.nan)
    metadata.index = pd.MultiIndex.from_tuples([sample_key(p, s) for p, s in zip(metadata.Project, metadata.SampleID)])
    assert len(set(keys)) == len(keys) and metadata.index.is_unique
    assert set(keys) == set(metadata.index), "Activity and metadata sample sets differ"
    metadata = metadata.loc[keys]
    assert np.isfinite(activity).all() and len(set(modules)) == len(modules)
    rng = np.random.default_rng(0)
    rows = []
    for field in FIELDS:
        labels = metadata[field]
        if field == "Time point (day)":
            labels = pd.qcut(pd.to_numeric(labels, errors="coerce"), q=4, duplicates="drop").astype(object)
        keep = labels.notna().to_numpy()
        values = activity[:, keep]
        labels = labels[keep].astype(str).to_numpy()
        observed = effect(values, labels)
        if observed is None:
            raise ValueError(f"Insufficient observations for {field}")
        omega, f_stat, p_anova, n, k = observed
        hits = np.zeros(len(modules), dtype=int)
        for _ in range(1000):
            permuted = labels.copy()
            rng.shuffle(permuted)
            hits += effect(values, permuted)[0] >= omega
        for j, module in enumerate(modules):
            rows.append(dict(iModulon=module, metadata=field, omega2=omega[j],
                f_stat=f_stat[j], p_anova=p_anova[j], p_perm=(hits[j] + 1) / 1001, n=n, k=k))
        print(f"{field}: n={n}, groups={k}", flush=True)
    table = pd.DataFrame(rows)
    table["q_value"] = stats.false_discovery_control(table.p_perm)
    table["FDR05"] = table.q_value < .05
    table.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
