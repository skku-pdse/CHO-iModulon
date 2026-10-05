"""Run the CHO TPM-to-iModulon workflow extracted from the author's ICA.py."""

import os

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"

import argparse
import hashlib
import importlib.metadata
import json
import logging
from pathlib import Path
import re
import time

import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler, normalize


def read_seeds(path):
    pairs = []
    for line in path.read_text().splitlines():
        match = re.fullmatch(r"\s*Run\s+(\d+)\s*:\s*Seed\s+(\d+)\s*", line)
        if match:
            pairs.append(tuple(map(int, match.groups())))
        elif line.strip():
            raise ValueError(f"Invalid seed line: {line}")
    pairs.sort()
    if not pairs or [i for i, _ in pairs] != list(range(1, len(pairs) + 1)):
        raise ValueError("Seed run indices must be consecutive, starting at 1.")
    return pairs


def preprocess_tpm(path):
    """Final, already sample-filtered TPM; standardize each sample across genes."""
    tpm = pd.read_csv(path, sep="\t", header=[0, 1], index_col=0)
    # Remove the discovery spreadsheet's unlabeled padding (including a stray
    # zero at the intersection of its final unlabeled row and column).
    unnamed = np.array([all(str(v).startswith("Unnamed:") for v in key)
                        for key in tpm.columns])
    no_gene = tpm.index.isna()
    if no_gene.any():
        if tpm.loc[no_gene, ~unnamed].notna().any().any():
            raise ValueError("Unlabeled gene row contains sample measurements.")
        tpm = tpm.loc[~no_gene]
    if unnamed.any() and tpm.loc[:, unnamed].notna().any().any():
        raise ValueError("Unlabeled sample column contains gene measurements.")
    tpm = tpm.dropna(axis=0, how="all").dropna(axis=1, how="all")
    values = tpm.to_numpy(dtype=float)
    if (not np.isfinite(values).all() or (values < 0).any()
            or not tpm.index.is_unique or not tpm.columns.is_unique):
        raise ValueError("TPM must be finite, nonnegative, and have unique gene/sample keys.")
    tpm = tpm.loc[tpm.gt(1).any(axis=1)]
    if tpm.empty:
        raise ValueError("No genes pass TPM > 1.")
    return tpm, StandardScaler().fit_transform(np.log2(tpm + 1))


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_matrices(tpm, transformed, runs, directory, prefix, components):
    from picard import Picard

    for run, seed in runs:
        s_path = directory / f"{prefix}_run{run:03d}_weights.csv"
        a_path = directory / f"{prefix}_run{run:03d}_mixing.csv"
        if s_path.exists() and a_path.exists():
            logging.info("Run %s: using completed matrices", run)
            continue
        if s_path.exists() or a_path.exists():
            raise RuntimeError(f"Incomplete run {run}; use a new output directory.")
        started = time.perf_counter()
        np.random.seed(seed)
        model = Picard(n_components=components, tol=1e-7,
                       max_iter=int(1e10), random_state=seed)
        weights = model.fit_transform(transformed)
        pd.DataFrame(weights, index=tpm.index).to_csv(s_path)
        pd.DataFrame(model.mixing_, index=tpm.columns).to_csv(a_path)
        logging.info("Run %s seed %s completed in %.2f s", run, seed,
                     time.perf_counter() - started)


def build_distances(runs, directory, distances, prefix):
    for i, _ in runs:
        left = None
        for j, _ in runs:
            if j < i:
                continue
            output = distances / f"{prefix}_similarity_{i:03d}_{j:03d}.npz"
            if output.exists():
                continue
            if left is None:
                left = normalize(pd.read_csv(
                    directory / f"{prefix}_run{i:03d}_weights.csv", index_col=0), axis=0)
            right = left if i == j else normalize(pd.read_csv(
                directory / f"{prefix}_run{j:03d}_weights.csv", index_col=0), axis=0)
            similarity = np.abs(left.T @ right)
            similarity[similarity < 0.5] = 0
            sparse.save_npz(output, sparse.coo_matrix(np.clip(similarity, 0, 1)))
        logging.info("Similarity blocks for run %s completed", i)


def aggregate_clusters(weights, mixing, labels):
    """Preserve original sign alignment, mean M, and signed-RMS A aggregation."""
    n_clusters = int(labels.max()) + 1
    if n_clusters == 0:
        raise ValueError("No stable clusters found.")
    m = pd.DataFrame(index=weights[0].index)
    a = pd.DataFrame(index=mixing[0].index)
    summaries = []
    flat_weights = pd.concat(weights, axis=1)
    flat_mixing = pd.concat(mixing, axis=1)
    for label in range(n_clusters):
        indices = np.flatnonzero(labels == label)
        s = flat_weights.iloc[:, indices].to_numpy(copy=True)
        activity = flat_mixing.iloc[:, indices].to_numpy(copy=True)
        if abs(s[:, 0].min()) > s[:, 0].max():
            s[:, 0] *= -1
            activity[:, 0] *= -1
        for j in range(1, s.shape[1]):
            if not stats.pearsonr(s[:, j], s[:, 0])[0] > 0:
                s[:, j] *= -1
                activity[:, j] *= -1
        m[label] = s.mean(axis=1)
        a[label] = np.sign(activity.mean(axis=1)) * np.sqrt((activity ** 2).mean(axis=1))
        summaries.append({"S_mean_std": s.std(axis=0).mean(),
                          "A_mean_std": activity.std(axis=0).mean(), "count": len(indices)})
    return m, a.T, pd.DataFrame(summaries)


def cluster_matrices(runs, directory, distances, prefix, output):
    blocks = []
    for i, _ in runs:
        row = []
        for j, _ in runs:
            low, high = sorted((i, j))
            block = sparse.load_npz(distances / f"{prefix}_similarity_{low:03d}_{high:03d}.npz")
            row.append(block if i <= j else block.T)
        blocks.append(row)
    similarity = sparse.bmat(blocks, format="csr")
    distance = sparse.csr_matrix((1 - similarity.data, similarity.indices,
                                  similarity.indptr), shape=similarity.shape)
    labels = DBSCAN(eps=0.1, min_samples=int(round(0.5 * (len(runs) + 1)) + 1),
                    metric="precomputed").fit_predict(distance)
    weights, mixing = [], []
    for run, _ in runs:
        weights.append(pd.read_csv(directory / f"{prefix}_run{run:03d}_weights.csv", index_col=0))
        mixing.append(pd.read_csv(directory / f"{prefix}_run{run:03d}_mixing.csv", index_col=[0, 1]))
        if not weights[-1].index.equals(weights[0].index) or not mixing[-1].index.equals(mixing[0].index):
            raise ValueError("Run matrices have inconsistent gene or sample order.")
    m, a, summary = aggregate_clusters(weights, mixing, labels)
    m.to_csv(output / f"{prefix}_ica_gene_weights.csv")
    a.to_csv(output / f"{prefix}_ica_activities.csv")
    summary.to_csv(output / f"{prefix}_ica_cluster_statistics.csv")
    logging.info("Saved %s iModulons", m.shape[1])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tpm", type=Path, required=True)
    parser.add_argument("--dataset", choices=["discovery", "expanded"], required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seeds", type=Path, default=Path(__file__).with_name("ica_seeds_256.txt"))
    parser.add_argument("--components", type=int, default=200)
    parser.add_argument("--stage", choices=["all", "ica", "distances", "cluster"], default="all")
    args = parser.parse_args()
    if args.components < 1:
        parser.error("--components must be positive")
    runs = read_seeds(args.seeds)
    versions = {name: importlib.metadata.version(name) for name in
                ("numpy", "pandas", "scipy", "scikit-learn", "python-picard")}
    config = {"tpm_sha256": sha256(args.tpm), "runs": runs, "components": args.components,
              "dataset": args.dataset, "versions": versions}
    config = json.loads(json.dumps(config))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = args.output_dir / f"{args.dataset}_ica_configuration.json"
    if manifest.exists():
        if json.loads(manifest.read_text()) != config:
            raise ValueError("Output directory belongs to a different input/configuration/environment.")
    elif any(args.output_dir.iterdir()):
        raise ValueError("Use an empty output directory for a new run.")
    else:
        manifest.write_text(json.dumps(config, indent=2) + "\n")
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", handlers=[
        logging.StreamHandler(), logging.FileHandler(args.output_dir / f"{args.dataset}_ica_execution.log")])
    matrices, distances = args.output_dir / "matrices", args.output_dir / "similarities"
    matrices.mkdir(exist_ok=True)
    distances.mkdir(exist_ok=True)
    logging.info("Stage: %s; input: %s", args.stage, args.tpm)
    if args.stage in ("all", "ica"):
        tpm, transformed = preprocess_tpm(args.tpm)
        if args.components > min(tpm.shape):
            raise ValueError("Components exceed input dimensions.")
        logging.info("Input after filtering: %s genes, %s samples", *tpm.shape)
        run_matrices(tpm, transformed, runs, matrices, args.dataset, args.components)
    if args.stage in ("all", "distances"):
        build_distances(runs, matrices, distances, args.dataset)
    if args.stage in ("all", "cluster"):
        cluster_matrices(runs, matrices, distances, args.dataset, args.output_dir)


if __name__ == "__main__":
    main()
