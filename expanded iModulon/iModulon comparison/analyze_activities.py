"""Current Results 3.5: reuse legacy statistical helpers; explicit keyed joins."""

from pathlib import Path
import json
import hashlib
import re
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image
from comparison_paths import DATA, RESULTS, read_weights
from statistics_helpers import (
    pearson as corr,
    project_center as center,
    exact_two_sided_permutation_p as perm,
)



BLUE = "#0072B2"
ORANGE = "#E69F00"
GREEN = "#009E73"
PURPLE = "#6F4AA8"
GRAY = "#BDBDBD"


def save_table(activity_comparison, n):
    activity_comparison.to_csv(RESULTS / n, index=False, encoding="utf-8-sig")


def normalize_sample_id(x):
    """Normalize the verified numeric suffix convention; uniqueness is checked later."""
    s = str(x).strip()
    if s.endswith(".0") and s[:-2].isdigit():
        s = s[:-2]
    return re.sub("\\.\\d+$", "", s)


def normalize_project(x):
    """Resolve the two verified aliases for the NIST project."""
    return {
        "R17704": "NISTCHO_Characterisation",
        "Project_BOKU": "NISTCHO_Characterisation",
    }.get(str(x).strip(), str(x).strip())


def sample_key(p, s):
    """Construct a composite sample key including the confirmed ZeLa ID correction."""
    p = normalize_project(p)
    s = normalize_sample_id(s)
    if p == "2020_02_12_ZeN" and s == "D2-Z-ZeLa-2_S17":
        s = "D3-Z-ZeLa-2_S17"
    return (p, s)


def read_activity(raw, label):
    """Read a two-header activity matrix and reject duplicate keys or nonfinite values."""
    ids = pd.to_numeric(raw.iloc[2:, 0], errors="coerce")
    good = ids.notna()
    assert raw.iloc[2:, 1:].loc[~good].isna().all().all(), (
        label + " unexpected nonnumeric rows"
    )
    values = raw.iloc[2:, 1:].loc[good].to_numpy(float)
    keys = [sample_key(p, s) for p, s in zip(raw.iloc[0, 1:], raw.iloc[1, 1:])]
    assert (
        len(set(keys)) == len(keys)
        and ids[good].is_unique
        and np.isfinite(values).all()
    ), label
    return pd.DataFrame(
        values,
        index=ids[good].astype(int),
        columns=pd.MultiIndex.from_tuples(keys, names=["Project", "SampleID"]),
    )


# Validate all inputs before calculating correlations.
print("Loading activity matrices", flush=True)
discovery_activity = read_activity(
    pd.read_csv(DATA / "discovery_A.csv.gz", header=None), "discovery"
)
expanded_activity = read_activity(
    pd.read_csv(DATA / "expanded_A.csv.gz", header=None), "expanded"
)
projected_activity = read_activity(
    pd.read_csv(DATA / "projected_A.csv.gz", header=None), "projected"
)
assert (
    discovery_activity.shape == (105, 778)
    and expanded_activity.shape == (148, 1030)
    and (projected_activity.shape == (105, 255))
    and (set(projected_activity.index) == set(discovery_activity.index))
)
assert set(discovery_activity.columns) <= set(expanded_activity.columns)
external = [
    k for k in expanded_activity.columns if k not in set(discovery_activity.columns)
]
assert len(external) == 252
vm = pd.read_csv(DATA / "projected_metadata.csv")
vm.columns = vm.columns.str.strip()
vk = [sample_key(p, s) for p, s in zip(vm.Project, vm.SampleID)]
assert len(set(vk)) == len(vk) and set(vk) == set(projected_activity.columns)
tech = [k for k in projected_activity.columns if k[0] == "PRJNA214684"]
assert len(tech) == 3 and set(tech) <= set(discovery_activity.columns)
assert set(projected_activity.columns) - set(tech) == set(external)
assert set((k[1] for k in vk if k[0] == "NISTCHO_Characterisation")) == set(
    (k[1] for k in expanded_activity.columns if k[0] == "NISTCHO_Characterisation")
)
projected_activity = projected_activity.loc[:, external]
shared_expanded_activity = expanded_activity.loc[:, discovery_activity.columns]
external_expanded_activity = expanded_activity.loc[:, external]
metadata = pd.read_csv(DATA / "metadata.csv")
metadata.columns = metadata.columns.str.strip()
mk = [sample_key(p, s) for p, s in zip(metadata.Project, metadata.SampleID)]
assert len(set(mk)) == len(mk) and set(mk) == set(expanded_activity.columns)
metadata.index = pd.MultiIndex.from_tuples(mk)
metadata = metadata.loc[expanded_activity.columns].copy()
save_table(
    pd.DataFrame(
        [
            {
                "cohort": (
                    "shared" if k in set(discovery_activity.columns) else "external"
                ),
                "Project": k[0],
                "SampleID": k[1],
            }
            for k in expanded_activity.columns
        ]
    ),
    "sample_alignment.csv",
)
best_matches = pd.read_csv(RESULTS / "discovery_best_matches.csv").set_index(
    "discovery_iModulon"
)
membership_pairs = pd.read_csv(RESULTS / "all_pairwise_membership.csv")
expanded_names = pd.read_csv(DATA / "expanded_names.csv").set_index("iModulon")
discovery_names = pd.read_csv(DATA / "discovery_names.csv").set_index("iModulon")
discovery_members = pd.read_csv(DATA / "discovery_members.csv")
expanded_members = pd.read_csv(DATA / "expanded_members.csv")
discovery_weights = read_weights("discovery")
expanded_weights = read_weights("expanded")
discovery_weights.columns = discovery_weights.columns.astype(int)
expanded_weights.columns = expanded_weights.columns.astype(int)
assert (
    discovery_weights.index.is_unique
    and expanded_weights.index.is_unique
    and (set(discovery_weights.index) == set(expanded_weights.index))
)
expanded_weights = expanded_weights.loc[discovery_weights.index]
# Compare activities only after keyed alignment and ICA sign correction.
shared_projects = np.array(discovery_activity.columns.get_level_values(0))
rows = []
projectrows = []
for d, b in best_matches.iterrows():
    e = int(b.expanded_iModulon)
    sign = np.sign(b.pearson_r)
    x = discovery_activity.loc[d].to_numpy()
    y = shared_expanded_activity.loc[e].to_numpy() * sign
    m = membership_pairs[
        membership_pairs.discovery_iModulon.eq(d)
        & membership_pairs.expanded_iModulon.eq(e)
    ].iloc[0]
    rows.append(
        dict(
            discovery_iModulon=d,
            expanded_iModulon=e,
            discovery_name=discovery_names.loc[d, "Name"],
            expanded_name=expanded_names.loc[e, "Name"],
            weight_r=b.pearson_r,
            weight_abs_r=b.abs_pearson_r,
            weight_RBM=bool(b.mutual_nearest_neighbor),
            membership_RBM=bool(m.reciprocal_unique),
            membership_J=m.membership_jaccard,
            discovery_n=m.discovery_n,
            expanded_n=m.expanded_n,
            shared_genes=m.shared_genes,
            activity_r=corr(x, y),
            centered_activity_r=corr(
                center(x, shared_projects), center(y, shared_projects)
            ),
            n_shared=len(x),
        )
    )
    for p in pd.unique(shared_projects):
        mask = shared_projects == p
        if mask.sum() >= 4:
            projectrows.append(
                dict(
                    discovery_iModulon=d,
                    expanded_iModulon=e,
                    Project=p,
                    n=int(mask.sum()),
                    activity_r=corr(x[mask], y[mask]),
                )
            )
activity_comparison = pd.DataFrame(rows)
reciprocal_pairs = activity_comparison[activity_comparison.weight_RBM]
assert len(reciprocal_pairs) == 89
save_table(activity_comparison, "activity_and_membership_comparison.csv")
save_table(pd.DataFrame(projectrows), "within_project_activity_correlations.csv")
summary = {
    "discovery_modules": 105,
    "expanded_modules": 148,
    "shared_samples": 778,
    "external_samples": len(external),
    "shared_genes": len(discovery_weights),
    "weight_RBM": len(reciprocal_pairs),
    "membership_unique_RBM": int(membership_pairs.reciprocal_unique.sum()),
    "concordant_RBM": int(
        (membership_pairs.reciprocal_unique & membership_pairs.weight_reciprocal).sum()
    ),
}
for label, col in [
    ("gene_weight", "weight_abs_r"),
    ("shared_activity", "activity_r"),
    ("centered_activity", "centered_activity_r"),
]:
    summary[label] = {
        "median": float(reciprocal_pairs[col].median()),
        "Q1": float(reciprocal_pairs[col].quantile(0.25)),
        "Q3": float(reciprocal_pairs[col].quantile(0.75)),
    }
selected = [87, 54, 39, 52, 65]
# Projection values are precomputed; no basis fitting occurs here.
externalrows = []
sample_rows = []
for d in [52, 65, 39, 54]:
    e = int(best_matches.loc[d, "expanded_iModulon"])
    sign = np.sign(best_matches.loc[d, "pearson_r"])
    x = projected_activity.loc[d].to_numpy()
    y = external_expanded_activity.loc[e].to_numpy() * sign
    externalrows.append(
        dict(
            discovery_iModulon=d,
            expanded_iModulon=e,
            Name=discovery_names.loc[d, "Name"],
            n=len(x),
            pearson_r=corr(x, y),
        )
    )
    for k, xx, yy in zip(external, x, y):
        sample_rows.append(
            dict(
                discovery_iModulon=d,
                expanded_iModulon=e,
                Project=k[0],
                SampleID=k[1],
                projected_activity=xx,
                expanded_sign_aligned_activity=yy,
            )
        )
ext = pd.DataFrame(externalrows)
save_table(ext, "external_projection_correlations.csv")
save_table(pd.DataFrame(sample_rows), "external_projection_sample_values.csv")
save_table(
    activity_comparison[activity_comparison.discovery_iModulon.isin(selected)],
    "representative_modules.csv",
)
overlap = []
for d in selected:
    e = int(best_matches.loc[d, "expanded_iModulon"])
    ds = set(discovery_members[discovery_members.iModulon.eq(d)].Gene)
    es = set(expanded_members[expanded_members.iModulon.eq(e)].Gene)
    for gene in sorted(ds | es):
        overlap.append(
            dict(
                discovery_iModulon=d,
                expanded_iModulon=e,
                Gene=gene,
                status=(
                    "shared"
                    if gene in ds & es
                    else "discovery_only" if gene in ds else "expanded_only"
                ),
                discovery_loading=discovery_weights.loc[gene, d],
                expanded_sign_aligned_loading=expanded_weights.loc[gene, e]
                * np.sign(best_matches.loc[d, "pearson_r"]),
            )
        )
save_table(pd.DataFrame(overlap), "representative_member_gene_changes.csv")
cell = metadata.copy()
cell["Activity"] = expanded_activity.loc[86].to_numpy()
cell["Cell line"] = cell["Cell line"].fillna("Unknown").str.strip()
save_table(cell.reset_index(drop=True), "expanded_iM86_sample_activities.csv")
group = (
    cell.groupby(["Project", "Cell line"])
    .Activity.agg(["count", "mean", "std"])
    .reset_index()
)
save_table(group, "expanded_iM86_project_cellline_summary.csv")
# Enumerate every label assignment for the specified within-project contrast.
within = cell[
    cell.Project.eq("PRJNA1130622") & cell["Cell line"].isin(["CHO-DXB11", "CHO-K1"])
]
x = within[within["Cell line"].eq("CHO-DXB11")].Activity.to_numpy()
y = within[within["Cell line"].eq("CHO-K1")].Activity.to_numpy()
assert (len(x), len(y)) == (22, 6)
print("Exact permutation: 376,740 label assignments", flush=True)
pv, nperm = perm(x, y)
assert nperm == 376740
pr = {
    "expanded_iModulon": 86,
    "Name": expanded_names.loc[86, "Name"],
    "Project": "PRJNA1130622",
    "n_DXB11": len(x),
    "n_K1": len(y),
    "mean_DXB11": float(x.mean()),
    "mean_K1": float(y.mean()),
    "difference_DXB11_minus_K1": float(x.mean() - y.mean()),
    "exact_two_sided_P": pv,
    "assignments": nperm,
}
save_table(pd.DataFrame([pr]), "DXB11_exact_permutation.csv")
summary["DXB11_test"] = pr
summary["external_correlations"] = externalrows


(RESULTS / "statistics_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
