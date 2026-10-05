"""Statistical helpers retained from the original analysis."""

import numpy as np
import pandas as pd
import itertools


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.std(x, ddof=1) == 0 or np.std(y, ddof=1) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def project_center(values: np.ndarray, projects: np.ndarray) -> np.ndarray:
    centered = np.empty_like(values, dtype=float)
    for project in pd.unique(projects):
        mask = projects == project
        centered[mask] = values[mask] - values[mask].mean()
    return centered


def exact_two_sided_permutation_p(x: np.ndarray, y: np.ndarray) -> tuple[float, int]:
    values = np.concatenate([x, y])
    observed = float(x.mean() - y.mean())
    n_x = len(x)
    extreme = 0
    total = 0
    for selected in itertools.combinations(range(len(values)), n_x):
        mask = np.zeros(len(values), dtype=bool)
        mask[list(selected)] = True
        difference = float(values[mask].mean() - values[~mask].mean())
        total += 1
        if abs(difference) >= abs(observed) - 1e-12:
            extreme += 1
    return (extreme / total, total)
