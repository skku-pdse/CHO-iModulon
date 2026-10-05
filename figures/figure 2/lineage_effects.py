"""Required generator for the pairwise lineage figure; v4 metadata used.

Adapted from lineage_effect_analysis.py. Original inputs remain unchanged.
"""

import matplotlib

matplotlib.use("Agg")
from figure2_paths import DATA, RESULTS
import re
from itertools import combinations
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

IMODULON_FILE = DATA / "discovery.xlsx"
METADATA_FILE = DATA / "metadata.xlsx"
OUTDIR = RESULTS / "lineage"
OUTDIR.mkdir(parents=True, exist_ok=True)
GROUP_ORDER = ["CHO-K1", "CHO-S", "CHO-DG44", "CHO-ZeLa"]
MIN_GROUP_N = 5
MIN_REST_N = 10
BOOT_N = 2000
BOOT_SEED = 20250425
BIMODAL_RATIO = 4.0


def normalize_sample_name(x: str) -> str:
    if pd.isna(x):
        return x
    x = str(x).strip()
    x = re.sub("\\s+", "", x)
    return x


def load_a_matrix(imodulon_file: str, sheet_name: str = "A") -> pd.DataFrame:
    A = pd.read_excel(imodulon_file, sheet_name=sheet_name, index_col=0, header=1)
    A.columns = [normalize_sample_name(c) for c in A.columns]
    A.index = [str(i) for i in A.index]
    return A


def load_metadata(metadata_file: str, sheet_name: str = "Main_v3") -> pd.DataFrame:
    meta = pd.read_excel(metadata_file, sheet_name=sheet_name)
    if "SampleID" not in meta.columns:
        raise ValueError(f"'SampleID' column not found. columns={list(meta.columns)}")
    meta = meta.copy()
    meta["Sample"] = meta["SampleID"].map(normalize_sample_name)
    # Current activity and metadata inputs already share the corrected D3 identifier.
    if meta["Sample"].duplicated().any():
        counts = {}
        new_names = []
        for s in meta["Sample"]:
            counts[s] = counts.get(s, 0) + 1
            new_names.append(s if counts[s] == 1 else f"{s}.{counts[s] - 1}")
        meta["Sample"] = new_names
    return meta


def merge_a_and_metadata(A: pd.DataFrame, meta: pd.DataFrame) -> pd.DataFrame:
    A_long = (
        A.T.reset_index()
        .rename(columns={"index": "Sample"})
        .melt(id_vars="Sample", var_name="iModulon", value_name="Activity")
    )
    assert meta.Sample.is_unique and A.columns.is_unique
    assert set(A.columns) == set(meta.Sample), "Unmatched sample identifiers"
    return A_long.merge(meta, on="Sample", how="inner", validate="many_to_one")


def cohens_d(x, y):
    """Cohen's d using pooled SD. Sign convention: x - y."""
    x = pd.Series(x).dropna().astype(float)
    y = pd.Series(y).dropna().astype(float)
    nx, ny = (len(x), len(y))
    if nx < 2 or ny < 2:
        return np.nan
    var_x = x.var(ddof=1)
    var_y = y.var(ddof=1)
    pooled_var = ((nx - 1) * var_x + (ny - 1) * var_y) / (nx + ny - 2)
    if pooled_var <= 0 or np.isnan(pooled_var):
        return np.nan
    return (x.mean() - y.mean()) / np.sqrt(pooled_var)


def hedges_g(x, y):
    """Hedges' g = small-sample-corrected Cohen's d. Use for n<20."""
    d = cohens_d(x, y)
    if np.isnan(d):
        return np.nan
    n = len(pd.Series(x).dropna()) + len(pd.Series(y).dropna())
    if n - 2 <= 0:
        return d
    j = 1 - 3 / (4 * (n - 2) - 1)
    return d * j


def auroc_two_groups(x, y):
    """
    AUROC of x vs y. Returns probability that a random x > random y.
    Equivalent to Mann-Whitney U / (n_x * n_y).
    """
    x = pd.Series(x).dropna().astype(float).values
    y = pd.Series(y).dropna().astype(float).values
    if len(x) < 2 or len(y) < 2:
        return np.nan
    try:
        u, _ = stats.mannwhitneyu(x, y, alternative="two-sided")
        return u / (len(x) * len(y))
    except ValueError:
        return np.nan


def compute_pairwise_effects(merged_df: pd.DataFrame, group_order: list = GROUP_ORDER):
    """
    For each iModulon, compute pairwise effect sizes for all 6 lineage pairs
    (or fewer if some are filtered). Sign convention: positive value means
    group_A > group_B.
    """
    df = merged_df[merged_df["Cell line"].isin(group_order)].copy()
    df = df.dropna(subset=["Activity"])
    counts = df.groupby("Cell line")["Sample"].nunique()
    valid_groups = [g for g in group_order if counts.get(g, 0) >= MIN_GROUP_N]
    pairs = list(combinations(valid_groups, 2))
    rows = []
    for imod, sub in df.groupby("iModulon"):
        for ga, gb in pairs:
            x = sub.loc[sub["Cell line"] == ga, "Activity"]
            y = sub.loc[sub["Cell line"] == gb, "Activity"]
            n_a, n_b = (len(x), len(y))
            if n_a < MIN_GROUP_N or n_b < MIN_GROUP_N:
                continue
            d = cohens_d(x, y)
            g = hedges_g(x, y)
            auc = auroc_two_groups(x, y)
            mean_diff = x.mean() - y.mean()
            rows.append(
                {
                    "iModulon": imod,
                    "group_A": ga,
                    "group_B": gb,
                    "n_A": int(n_a),
                    "n_B": int(n_b),
                    "mean_A": x.mean(),
                    "mean_B": y.mean(),
                    "sd_A": x.std(ddof=1),
                    "sd_B": y.std(ddof=1),
                    "mean_diff": mean_diff,
                    "cohens_d": d,
                    "hedges_g": g,
                    "auroc": auc,
                    "auroc_disc": np.nan if np.isnan(auc) else abs(auc - 0.5),
                }
            )
    return pd.DataFrame(rows)


def oneway_anova_with_etasq(merged_df: pd.DataFrame, group_order: list = GROUP_ORDER):
    """
    One-way ANOVA across cell lines (n >= MIN_GROUP_N) per iModulon.
    Reports F, p, eta_sq, omega_sq (less biased for small samples).
    """
    df = merged_df[merged_df["Cell line"].isin(group_order)].copy()
    df = df.dropna(subset=["Activity"])
    counts = df.groupby("Cell line")["Sample"].nunique()
    valid_groups = [g for g in group_order if counts.get(g, 0) >= MIN_GROUP_N]
    df = df[df["Cell line"].isin(valid_groups)].copy()
    rows = []
    for imod, sub in df.groupby("iModulon"):
        groups = [
            sub.loc[sub["Cell line"] == g, "Activity"].values for g in valid_groups
        ]
        ns = [len(v) for v in groups]
        if any((n < MIN_GROUP_N for n in ns)):
            continue
        f_stat, p_val = stats.f_oneway(*groups)
        all_vals = np.concatenate(groups)
        grand = all_vals.mean()
        ss_between = sum((n * (g.mean() - grand) ** 2 for g, n in zip(groups, ns)))
        ss_within = sum((((g - g.mean()) ** 2).sum() for g in groups))
        ss_total = ss_between + ss_within
        k = len(valid_groups)
        n_total = sum(ns)
        df_between = k - 1
        df_within = n_total - k
        eta_sq = ss_between / ss_total if ss_total > 0 else np.nan
        ms_within = ss_within / df_within if df_within > 0 else np.nan
        omega_sq = (
            (ss_between - df_between * ms_within) / (ss_total + ms_within)
            if ss_total + ms_within > 0
            else np.nan
        )
        rows.append(
            {
                "iModulon": imod,
                "anova_F": f_stat,
                "anova_p": p_val,
                "eta_sq": eta_sq,
                "omega_sq": omega_sq,
                "n_total": int(n_total),
                "k_groups": int(k),
            }
        )
    out = pd.DataFrame(rows)
    out["anova_p_bonf"] = (out["anova_p"] * len(out)).clip(upper=1.0)
    return out


def bootstrap_lineage_means(
    merged_df: pd.DataFrame,
    group_order: list = GROUP_ORDER,
    n_boot: int = BOOT_N,
    seed: int = BOOT_SEED,
):
    """
    For each iModulon x lineage, resample with replacement n_boot times.
    Report point estimate, 2.5/97.5 percentiles, bootstrap SE.
    Primarily defends small-n lineages (DG44 n=10, ZeLa n=29).
    """
    rng = np.random.default_rng(seed)
    df = merged_df[merged_df["Cell line"].isin(group_order)].copy()
    df = df.dropna(subset=["Activity"])
    rows = []
    for imod, sub in df.groupby("iModulon"):
        for g in group_order:
            vals = sub.loc[sub["Cell line"] == g, "Activity"].values
            n = len(vals)
            if n < MIN_GROUP_N:
                continue
            boots = np.empty(n_boot)
            for i in range(n_boot):
                idx = rng.integers(0, n, n)
                boots[i] = vals[idx].mean()
            mean_pt = vals.mean()
            ci_lo, ci_hi = np.percentile(boots, [2.5, 97.5])
            se = boots.std(ddof=1)
            rows.append(
                {
                    "iModulon": imod,
                    "group": g,
                    "n": int(n),
                    "mean": mean_pt,
                    "ci_lo": ci_lo,
                    "ci_hi": ci_hi,
                    "ci_width": ci_hi - ci_lo,
                    "boot_se": se,
                }
            )
    return pd.DataFrame(rows)


def within_lineage_variability(
    merged_df: pd.DataFrame, group_order: list = GROUP_ORDER
):
    """
    Per-lineage within-group SD per iModulon. Flag any lineage whose within-SD
    is >= BIMODAL_RATIO * min within-SD of other lineages in the same iModulon.
    Catches "K1 sub-population bimodality" patterns (iM25, iM88) that would be
    hidden in a mean heatmap.
    """
    df = merged_df[merged_df["Cell line"].isin(group_order)].copy()
    df = df.dropna(subset=["Activity"])
    rows = []
    for imod, sub in df.groupby("iModulon"):
        sds = {}
        ns = {}
        for g in group_order:
            vals = sub.loc[sub["Cell line"] == g, "Activity"]
            if len(vals) < MIN_GROUP_N:
                sds[g] = np.nan
                ns[g] = 0
            else:
                sds[g] = float(vals.std(ddof=1))
                ns[g] = len(vals)
        valid_sds = {g: s for g, s in sds.items() if not np.isnan(s)}
        if len(valid_sds) < 2:
            continue
        sd_min = min(valid_sds.values())
        sd_min_safe = sd_min if sd_min > 0 else 1e-09
        for g, s in sds.items():
            flag = not np.isnan(s) and s >= BIMODAL_RATIO * sd_min_safe and (s > 1.0)
            rows.append(
                {
                    "iModulon": imod,
                    "group": g,
                    "n": int(ns[g]),
                    "within_sd": s,
                    "min_within_sd_other_lineages": sd_min,
                    "bimodal_flag": bool(flag),
                }
            )
    return pd.DataFrame(rows)


def signed_dominant_lineage(
    pairwise_df: pd.DataFrame,
    mean_act_df: pd.DataFrame,
    group_order: list = GROUP_ORDER,
):
    """
    Identify the lineage whose mean is most distant from the centroid of the
    other lineages, and label direction ('gain' = higher than rest,
    'loss' = lower than rest). Avoids OvR misleading "dominant_group" labels.
    """
    rows = []
    for imod in mean_act_df.index:
        means = mean_act_df.loc[imod, group_order]
        if means.isna().any():
            continue
        best_g = None
        best_dist = np.nan
        best_dir = None
        for g in group_order:
            other = [m for og, m in means.items() if og != g]
            other_centroid = np.mean(other)
            dist = means[g] - other_centroid
            if best_g is None or abs(dist) > abs(best_dist):
                best_dist = dist
                best_g = g
                best_dir = "gain" if dist > 0 else "loss"
        p_sub = pairwise_df[
            (pairwise_df["iModulon"] == imod)
            & ((pairwise_df["group_A"] == best_g) | (pairwise_df["group_B"] == best_g))
        ]
        if len(p_sub):
            max_abs_d = p_sub["cohens_d"].abs().max()
            max_auroc_disc = p_sub["auroc_disc"].max()
        else:
            max_abs_d = np.nan
            max_auroc_disc = np.nan
        rows.append(
            {
                "iModulon": imod,
                "outlier_lineage": best_g,
                "direction": best_dir,
                "lineage_minus_others_centroid": best_dist,
                "max_abs_pairwise_d_for_outlier": max_abs_d,
                "max_auroc_disc_for_outlier": max_auroc_disc,
            }
        )
    return pd.DataFrame(rows)


def build_summary(
    pairwise_df: pd.DataFrame,
    anova_df: pd.DataFrame,
    bimodal_df: pd.DataFrame,
    signed_df: pd.DataFrame,
    mean_act_df: pd.DataFrame,
    group_order: list = GROUP_ORDER,
):
    """One row per iModulon. Suitable for ranking and direct heatmap input."""
    pw_max = (
        pairwise_df.groupby("iModulon")
        .agg(
            max_abs_d_pairwise=("cohens_d", lambda s: s.abs().max()),
            max_abs_g_pairwise=("hedges_g", lambda s: s.abs().max()),
            max_auroc_disc=("auroc_disc", "max"),
        )
        .reset_index()
    )

    def _max_pair(sub):
        sub = sub.dropna(subset=["cohens_d"])
        if sub.empty:
            return pd.Series(
                {"max_d_pair_A": np.nan, "max_d_pair_B": np.nan, "max_d_signed": np.nan}
            )
        i = sub["cohens_d"].abs().idxmax()
        return pd.Series(
            {
                "max_d_pair_A": sub.loc[i, "group_A"],
                "max_d_pair_B": sub.loc[i, "group_B"],
                "max_d_signed": sub.loc[i, "cohens_d"],
            }
        )

    pw_pair = (
        pairwise_df.groupby("iModulon", group_keys=False)
        .apply(_max_pair, include_groups=False)
        .reset_index()
    )
    bm_summary = (
        bimodal_df.groupby("iModulon")
        .agg(
            any_bimodal_flag=("bimodal_flag", "any"),
            n_bimodal_lineages=("bimodal_flag", "sum"),
            max_within_sd=("within_sd", "max"),
        )
        .reset_index()
    )
    out = anova_df.merge(pw_max, on="iModulon", how="left")
    out = out.merge(pw_pair, on="iModulon", how="left")
    out = out.merge(signed_df, on="iModulon", how="left")
    out = out.merge(bm_summary, on="iModulon", how="left")
    means_w = mean_act_df.copy()
    means_w.columns = [f"mean_{c}" for c in means_w.columns]
    means_w = means_w.reset_index().rename(
        columns={means_w.index.name or "index": "iModulon"}
    )
    if "iModulon" not in means_w.columns:
        means_w = means_w.rename(columns={means_w.columns[0]: "iModulon"})
    out = out.merge(means_w, on="iModulon", how="left")
    out = out.sort_values("eta_sq", ascending=False)
    return out


def build_mean_matrices(merged_df: pd.DataFrame, group_order: list = GROUP_ORDER):
    """Return (raw_means, zscore_means, within_sd_wide). Rows = iModulon."""
    df = merged_df[merged_df["Cell line"].isin(group_order)].copy()
    df = df.dropna(subset=["Activity"])
    means = (
        df.groupby(["iModulon", "Cell line"])["Activity"].mean().unstack("Cell line")
    )
    sds = (
        df.groupby(["iModulon", "Cell line"])["Activity"]
        .std(ddof=1)
        .unstack("Cell line")
    )
    means = means.reindex(columns=group_order)
    sds = sds.reindex(columns=group_order)
    means.index.name = "iModulon"
    sds.index.name = "iModulon"
    z = means.sub(means.mean(axis=1), axis=0).div(
        means.std(axis=1, ddof=1).replace(0, np.nan), axis=0
    )
    return (means, z, sds)


if __name__ == "__main__":
    print("[INFO] Loading data ...")
    A = load_a_matrix(IMODULON_FILE, sheet_name="A")
    meta = load_metadata(METADATA_FILE, sheet_name="Main_v3")
    merged = merge_a_and_metadata(A, meta)
    print(
        f"[INFO] Merged: n samples = {merged['Sample'].nunique()}, n iModulons = {merged['iModulon'].nunique()}"
    )
    print("[INFO] (1) Pairwise effect sizes ...")
    pairwise_df = compute_pairwise_effects(merged, GROUP_ORDER)
    pairwise_df.to_csv(OUTDIR / "lineage_pairwise_effect.csv", index=False)
    print(f"        -> {len(pairwise_df)} rows")
    print("[INFO] (2) ANOVA + eta_sq + omega_sq ...")
    anova_df = oneway_anova_with_etasq(merged, GROUP_ORDER)
    print("[INFO] (3) Bootstrap CI for lineage means ...")
    boot_df = bootstrap_lineage_means(
        merged, GROUP_ORDER, n_boot=BOOT_N, seed=BOOT_SEED
    )
    boot_df.to_csv(OUTDIR / "lineage_bootstrap_ci.csv", index=False)
    print(f"        -> {len(boot_df)} rows ({BOOT_N} bootstrap iterations)")
    print("[INFO] (4) Within-lineage SD + bimodality flag ...")
    bimodal_df = within_lineage_variability(merged, GROUP_ORDER)
    bimodal_df.to_csv(OUTDIR / "lineage_within_sd.csv", index=False)
    n_flag = bimodal_df["bimodal_flag"].sum()
    print(f"        -> {n_flag} lineage-iModulon cells flagged as bimodal")
    print("[INFO] (5) Mean activity matrices (raw + z-scored) ...")
    mean_act, z_act, sd_act = build_mean_matrices(merged, GROUP_ORDER)
    mean_act.to_csv(OUTDIR / "lineage_mean_activity.csv")
    z_act.to_csv(OUTDIR / "lineage_mean_zscore.csv")
    print("[INFO] (6) Sign-aware dominant lineage ...")
    signed_df = signed_dominant_lineage(pairwise_df, mean_act, GROUP_ORDER)
    print("[INFO] (7) Building summary table ...")
    summary_df = build_summary(
        pairwise_df, anova_df, bimodal_df, signed_df, mean_act, GROUP_ORDER
    )
    summary_df.to_csv(OUTDIR / "lineage_effect_summary.csv", index=False)
    print("\n=========================================================")
    print("Top 20 iModulons by eta_sq (cell line as factor):")
    cols = [
        "iModulon",
        "outlier_lineage",
        "direction",
        "eta_sq",
        "omega_sq",
        "anova_p_bonf",
        "max_abs_d_pairwise",
        "max_auroc_disc",
        "any_bimodal_flag",
        "max_d_pair_A",
        "max_d_pair_B",
    ]
    print(summary_df[cols].head(20).to_string(index=False))
    print("=========================================================\n")
    print(f"[DONE] All outputs saved to {OUTDIR}")
