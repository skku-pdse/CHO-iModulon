"""Adapt verified current tables to the unchanged original Figure 5 interfaces."""

from pathlib import Path
from functools import lru_cache
import pandas as pd
from figure5_paths import DATA, RESULTS, read_weights

S = RESULTS


def load_best_matches():
    return pd.read_csv(RESULTS / "discovery_best_matches.csv")


def load_activity():
    return pd.read_csv(S / "activity_and_membership_comparison.csv").rename(
        columns={
            "weight_r": "M_pearson_r_signed",
            "weight_abs_r": "M_abs_pearson_r",
            "weight_RBM": "mutual_nearest_neighbor",
            "activity_r": "activity_pearson_r_sign_aligned",
            "centered_activity_r": "within_project_centered_pearson_r",
        }
    )


@lru_cache(None)
def load_matched_data():
    return {
        "activity": load_activity(),
        "discovery_m": read_weights("discovery"),
        "expanded_m": read_weights("expanded"),
        "discovery_members": pd.read_csv(DATA / "discovery_members.csv"),
        "expanded_members": pd.read_csv(DATA / "expanded_members.csv"),
    }


@lru_cache(None)
def load_lineage_data():
    g = pd.read_csv(DATA / "expanded_members.csv")
    a = pd.read_csv(S / "expanded_iM86_sample_activities.csv").rename(
        columns={"Activity": "iM86_activity"}
    )
    return {
        "im86_members": g[g.iModulon.eq(86)],
        "im86_samples": a,
        "permutation_p": float(
            pd.read_csv(S / "DXB11_exact_permutation.csv").exact_two_sided_P.iloc[0]
        ),
    }
