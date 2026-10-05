"""Chromosome enrichment with a membership-union background and global BH FDR."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "expanded_chromosome_enrichment_results.csv")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(HERE.parent / "data/expanded_gene_membership.csv")
    assert not genes[["iModulon", "Gene", "Chromosome"]].isna().any().any()
    assert not genes.duplicated(["iModulon", "Gene"]).any()
    genes["Chromosome"] = genes.Chromosome.astype(str).replace({"0": "U", "Z": "U"})
    order = [str(i) for i in range(1, 11)] + ["X", "U"]
    assert set(genes.Chromosome) <= set(order)
    assert genes.groupby("Gene").Chromosome.nunique().max() == 1
    background = genes.drop_duplicates("Gene")
    total = len(background)
    counts = background.Chromosome.value_counts()
    rows = []
    for module, group in genes.groupby("iModulon", sort=True):
        size = len(group)
        if size < 5:
            rows.append(dict(iModulon=module, n_iModulon=size, status="skipped_fewer_than_5_genes"))
            continue
        for chromosome in order:
            successes = int(counts.get(chromosome, 0))
            if successes == 0:
                continue
            observed = int(group.Chromosome.eq(chromosome).sum())
            expected = size * successes / total
            rows.append(dict(iModulon=module, Chromosome=chromosome,
                n_iModulon=size, K_chromosome=successes, N_background=total,
                observed=observed, expected=expected,
                log2_fold=np.log2(observed / expected) if observed else np.nan,
                p_value=stats.hypergeom.sf(observed - 1, total, successes, size), status="tested"))
    table = pd.DataFrame(rows)
    tested = table.status.eq("tested")
    table.loc[tested, "p_adj_BH"] = stats.false_discovery_control(table.loc[tested, "p_value"])
    table["positive_enrichment_FDR05"] = table.p_adj_BH.lt(.05) & table.log2_fold.gt(0)
    table.to_csv(args.output, index=False)
    print(f"{tested.sum()} tests; background={total}; output=expanded_chromosome_enrichment_results.csv")


if __name__ == "__main__":
    main()
