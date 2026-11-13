"""End-to-end contract for the opt-in Gaussian best-flat-cut reference sibling p-value.

Same setting as ``65_test_split_count_rule_out_runner.py``: 30 coordinates and the registered
family-diffusion branch-length settings.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist

from benchmarks.shared.runners.dispatch import _tbs_branch_length_optimization_kwargs
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from benchmarks.shared.runners.tbs_runner import run_tbs_on_distance
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.best_flat_cut_reference import (
    MIN_CELLS,
)
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns

P = 30


def _fit(points: np.ndarray):
    params = dict(METHOD_SPECS["tbs_family_diffusion"].param_grid[0])
    data = pd.DataFrame(points, columns=[f"PC{i + 1}" for i in range(P)])
    return run_tbs_on_distance(
        data,
        pdist(points),
        0.01,
        edge_alpha=0.001,
        tree_linkage_method="ward",
        feature_space=continuous_feature_space_from_columns(tuple(data.columns)),
        branch_length_data_df=data,
        edge_branch_length_variance_policy=params["edge_branch_length_variance_policy"],
        **_tbs_branch_length_optimization_kwargs(params),
        sibling_best_flat_cut_reference=True,
    )


def test_reference_keeps_separated_groups_and_one_population() -> None:
    rng = np.random.default_rng(3)
    shift = np.r_[8.0, np.zeros(P - 1)]
    two_groups = np.vstack(
        [rng.standard_normal((150, P)) + shift, rng.standard_normal((150, P)) - shift]
    )
    result = _fit(two_groups)
    assert result.status == "ok"
    assert result.found_clusters == 2
    assert result.extra["sibling_best_flat_cut_reference"] is True
    assert result.extra["sibling_split_count_rule_out"] is False
    # The opt-in lets a peel below the reference's smallest tested size pass the edge gate.
    assert result.extra["gate_bundle"].metadata.config.small_peel_leaf_limit == MIN_CELLS
    methods = set(result.extra["annotations"]["Sibling_Test_Method"].dropna().astype(str))
    assert "best_flat_cut_reference:tested" in methods

    one_population = rng.standard_normal((300, P))
    assert _fit(one_population).found_clusters == 1
