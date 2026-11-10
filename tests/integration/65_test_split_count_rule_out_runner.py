"""End-to-end contract for the opt-in split-count rule-out sibling p-value.

Uses 30 coordinates and the registered family-diffusion branch-length settings, the setting the rule-out was
measured in; with fewer coordinates small subtrees keep invariant directions and noise splits can pass.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist

from benchmarks.shared.runners.dispatch import _tbs_branch_length_optimization_kwargs
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from benchmarks.shared.runners.tbs_runner import run_tbs_on_distance
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns

P = 30


def _fit(points: np.ndarray, *, rule_out: bool):
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
        sibling_split_count_rule_out=rule_out,
    )


def test_rule_out_keeps_separated_groups_and_one_population() -> None:
    rng = np.random.default_rng(3)
    shift = np.r_[8.0, np.zeros(P - 1)]
    two_groups = np.vstack(
        [rng.standard_normal((150, P)) + shift, rng.standard_normal((150, P)) - shift]
    )
    result = _fit(two_groups, rule_out=True)
    assert result.status == "ok"
    assert result.found_clusters == 2
    assert result.extra["sibling_split_count_rule_out"] is True
    methods = set(result.extra["annotations"]["Sibling_Test_Method"].dropna().astype(str))
    assert "split_count_rule_out:kept" in methods

    one_population = rng.standard_normal((300, P))
    assert _fit(one_population, rule_out=True).found_clusters == 1
    assert _fit(one_population, rule_out=False).found_clusters > 1  # today's test splits this noise
