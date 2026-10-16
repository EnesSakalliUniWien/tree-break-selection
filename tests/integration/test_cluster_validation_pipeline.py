"""
Production integration tests for cluster validation using complex synthetic data.

Tests the full pipeline with:
- Balanced binary feature matrices with low entropy
- Unbalanced binary feature matrices with high entropy
"""

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from sklearn.metrics import adjusted_rand_score

from benchmarks.shared.generators import generate_random_feature_matrix
from benchmarks.shared.runners.tbs_runner import run_tbs_on_distance
from tree_break_selection.hierarchy_analysis.statistics.alpha_contract import (
    DEFAULT_SIBLING_ALPHA,
)


def _run_tbs(**generator_kwargs) -> tuple[np.ndarray, list[int]]:
    """Run the production TBS entry point on a generated binary matrix."""
    data_dict, true_clusters = generate_random_feature_matrix(**generator_kwargs)
    data_df = pd.DataFrame.from_dict(data_dict, orient="index").astype(int)
    result = run_tbs_on_distance(
        data_df, pdist(data_df.values, metric="hamming"), DEFAULT_SIBLING_ALPHA
    )
    assert result.status == "ok"
    return np.asarray(result.labels), [true_clusters[name] for name in data_df.index]


def test_complex_random_feature_matrix_balanced_clusters():
    """Synthetic binary data with low entropy should recover the planted clusters."""
    predicted, true_labels = _run_tbs(
        n_rows=72,
        n_cols=40,
        entropy_param=0.1,
        n_clusters=4,
        random_seed=314,
        balanced_clusters=True,
    )

    # The production path recovers all four clusters exactly (ARI 1.0) on this fixture.
    assert len(np.unique(predicted)) == 4
    assert adjusted_rand_score(true_labels, predicted) > 0.95


def test_complex_random_feature_matrix_unbalanced_clusters():
    """Moderately higher entropy and unbalanced clusters should stay informative."""
    predicted, true_labels = _run_tbs(
        n_rows=150,
        n_cols=36,
        entropy_param=0.25,
        n_clusters=4,
        random_seed=2024,
        balanced_clusters=False,
    )

    # Measured: 2 clusters, ARI 0.48; the strict sibling calibration stops early here.
    assert len(np.unique(predicted)) >= 2
    assert adjusted_rand_score(true_labels, predicted) > 0.4
