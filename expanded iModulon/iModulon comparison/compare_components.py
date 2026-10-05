"""Compare ICA components using full weights and threshold-defined membership.

Top-100 overlap is retained as a historical reference, not a preservation cutoff.
Membership ties use exact fractions rather than rounded floating-point scores.
"""

from fractions import Fraction
import json
import numpy as np
import pandas as pd
from comparison_paths import DATA, RESULTS, read_weights


def main():
    discovery = read_weights("discovery")
    expanded = read_weights("expanded")
    for matrix in [discovery, expanded]:
        assert matrix.index.is_unique and matrix.columns.is_unique
        assert np.isfinite(matrix.to_numpy()).all()
        assert (matrix.std() > 0).all()
    shared = discovery.index.intersection(expanded.index)
    discovery = discovery.loc[shared]
    expanded = expanded.loc[shared]
    assert len(shared) == 19812

    # Center and normalize each component across genes to obtain Pearson r.
    def unit_vectors(matrix):
        values = matrix.to_numpy()
        values = values - values.mean(axis=0)
        return values / np.linalg.norm(values, axis=0)

    correlations = pd.DataFrame(
        unit_vectors(discovery).T @ unit_vectors(expanded),
        index=discovery.columns,
        columns=expanded.columns,
    )
    absolute = correlations.abs()
    forward = absolute.idxmax(axis=1)
    reverse = absolute.idxmax(axis=0)
    assert (
        np.isclose(
            absolute, absolute.max(axis=1).to_numpy()[:, None], atol=1e-12, rtol=0
        ).sum(axis=1)
        == 1
    ).all()
    assert (
        np.isclose(
            absolute, absolute.max(axis=0).to_numpy()[None, :], atol=1e-12, rtol=0
        ).sum(axis=0)
        == 1
    ).all()
    best = pd.DataFrame(
        [
            dict(
                discovery_iModulon=i,
                expanded_iModulon=j,
                pearson_r=correlations.loc[i, j],
                abs_pearson_r=absolute.loc[i, j],
                mutual_nearest_neighbor=bool(reverse[j] == i),
            )
            for i, j in forward.items()
        ]
    )
    best.to_csv(RESULTS / "discovery_best_matches.csv", index=False)
    pd.DataFrame(
        [
            dict(
                expanded_iModulon=j,
                discovery_iModulon=i,
                maximum_abs_pearson_r=absolute.loc[i, j],
            )
            for j, i in reverse.items()
        ]
    ).to_csv(RESULTS / "expanded_best_matches.csv", index=False)
    correlations.to_csv(RESULTS / "pairwise_gene_weight_correlations.csv")

    def member_sets(cohort):
        table = pd.read_csv(DATA / f"{cohort}_members.csv")
        assert not table[["iModulon", "Gene"]].isna().any().any()
        assert not table.duplicated(["iModulon", "Gene"]).any()
        return {i: set(group.Gene) for i, group in table.groupby("iModulon")}

    ds, es = member_sets("discovery"), member_sets("expanded")
    assert set(ds) == set(discovery.columns) and set(es) == set(expanded.columns)
    scores = {
        (i, j): Fraction(len(a & b), len(a | b))
        for i, a in ds.items()
        for j, b in es.items()
    }
    dmax = {i: max(scores[i, j] for j in es) for i in ds}
    emax = {j: max(scores[i, j] for i in ds) for j in es}
    dties = {i: sum(scores[i, j] == dmax[i] for j in es) for i in ds}
    eties = {j: sum(scores[i, j] == emax[j] for i in ds) for j in es}
    rows = []
    for (i, j), score in scores.items():
        mutual = score > 0 and score == dmax[i] and score == emax[j]
        rows.append(
            dict(
                discovery_iModulon=i,
                expanded_iModulon=j,
                discovery_n=len(ds[i]),
                expanded_n=len(es[j]),
                shared_genes=len(ds[i] & es[j]),
                union_genes=len(ds[i] | es[j]),
                membership_jaccard=float(score),
                membership_distance=1 - float(score),
                discovery_retained_fraction=len(ds[i] & es[j]) / len(ds[i]),
                expanded_covered_fraction=len(ds[i] & es[j]) / len(es[j]),
                reciprocal_including_ties=mutual,
                reciprocal_unique=mutual and dties[i] == eties[j] == 1,
                discovery_best_tie_count=dties[i],
                expanded_best_tie_count=eties[j],
                weight_reciprocal=bool(forward[i] == j and reverse[j] == i),
            )
        )
    pairs = pd.DataFrame(rows)
    pairs.to_csv(RESULTS / "all_pairwise_membership.csv", index=False)
    pairs[pairs.reciprocal_unique].to_csv(
        RESULTS / "membership_reciprocal_matches.csv", index=False
    )
    pairs[pairs.reciprocal_including_ties & ~pairs.reciprocal_unique].to_csv(
        RESULTS / "membership_tied_matches.csv", index=False
    )

    # Preserve the earlier top-100 comparison as a clearly labelled reference.
    top_d = {i: set(discovery[i].abs().nlargest(100).index) for i in discovery}
    top_e = {j: set(expanded[j].abs().nlargest(100).index) for j in expanded}
    historical = best.copy()
    historical["top_100_jaccard"] = [
        len(top_d[i] & top_e[j]) / len(top_d[i] | top_e[j])
        for i, j in zip(best.discovery_iModulon, best.expanded_iModulon)
    ]
    historical.to_csv(RESULTS / "historical_top100_reference.csv", index=False)
    summary = dict(
        weight_reciprocal_matches=int(best.mutual_nearest_neighbor.sum()),
        membership_reciprocal_matches=int(pairs.reciprocal_unique.sum()),
        concordant_pairs=int((pairs.reciprocal_unique & pairs.weight_reciprocal).sum()),
    )
    (RESULTS / "matching_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(summary)


if __name__ == "__main__":
    main()
