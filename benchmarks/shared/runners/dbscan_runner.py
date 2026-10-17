"""DBSCAN runner for precomputed distance matrices."""

from __future__ import annotations

import numpy as np
from sklearn.cluster import DBSCAN

from benchmarks.shared.types.method_run_result import MethodRunResult
from benchmarks.shared.util.core import estimate_dbscan_eps, normalize_labels
from benchmarks.shared.util.decomposition import ok_result_from_labels


def _run_dbscan_method(
    distance_matrix: np.ndarray,
    params: dict[str, object],
) -> MethodRunResult:
    """Run DBSCAN on a precomputed distance matrix."""

    n_samples = distance_matrix.shape[0]
    if n_samples <= 1:
        labels = np.zeros(n_samples, dtype=int)
        return ok_result_from_labels(labels, range(n_samples))

    min_samples = int(params["min_samples"])
    eps = params["eps"]
    if eps == "median_k_distance":
        eps = estimate_dbscan_eps(distance_matrix, min_samples)
    model = DBSCAN(metric="precomputed", eps=float(eps), min_samples=min_samples)
    labels = normalize_labels(model.fit_predict(distance_matrix))
    return ok_result_from_labels(labels, range(n_samples))
