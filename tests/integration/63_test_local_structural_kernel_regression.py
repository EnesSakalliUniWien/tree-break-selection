from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import pdist

from benchmarks.shared.cases import get_default_test_cases
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from benchmarks.shared.runners.tbs_runner import run_tbs_on_distance
from benchmarks.shared.util.case_inputs import prepare_case_inputs
from tree_break_selection.hierarchy_analysis.statistics.alpha_contract import (
    DEFAULT_SIBLING_ALPHA,
)
from tree_break_selection.hierarchy_analysis.statistics.branch_length_utils import (
    EDGE_BRANCH_LENGTH_VARIANCE_POLICY_NORMALIZED,
)


def _case_tbs_inputs(case_name: str):
    """Return a default case's data, TBS tree distance and feature space."""
    case = next(case for case in get_default_test_cases() if case["name"] == case_name)
    inputs = prepare_case_inputs(case, ["tbs"])
    metric = METHOD_SPECS["tbs"].param_grid[0]["tree_distance_metric"]
    return (
        inputs.data,
        pdist(inputs.data.values, metric=metric),
        inputs.metadata.get("feature_space"),
    )


@pytest.mark.slow
def test_null_case_returns_one_cluster() -> None:
    # A same-leaf edge basis rejected every null edge, which left no null-like sibling
    # record to calibrate the empirical null and made this case unsupported.
    data, distance, _feature_space = _case_tbs_inputs("gauss_null_large")

    result = run_tbs_on_distance(data, distance, DEFAULT_SIBLING_ALPHA)

    assert result.status == "ok"
    assert result.unsupported_reason is None
    assert result.found_clusters == 1


@pytest.mark.slow
def test_default_sibling_calibration_recovers_gauss_clear_small() -> None:
    data, distance, _feature_space = _case_tbs_inputs("gauss_clear_small")

    result = run_tbs_on_distance(
        data,
        distance,
        DEFAULT_SIBLING_ALPHA,
        edge_branch_length_variance_policy=EDGE_BRANCH_LENGTH_VARIANCE_POLICY_NORMALIZED,
        allow_linkage_ultrametric_branch_time=True,
        trace_level="full",
    )

    assert result.found_clusters == 3

    stage_timings = result.extra["stage_timings"]
    for key in (
        "tree_build_sec",
        "populate_divergences_sec",
        "edge_gate_sec",
        "edge_gate_contrast_covariance_sec",
        "edge_gate_projection_sec",
        "edge_gate_wald_statistic_sec",
        "edge_gate_tree_bh_sec",
        "spectral_context_sec",
        "tangent_whitening_sec",
        "eigensolve_sec",
        "pca_projection_sec",
        "sibling_gate_sec",
        "sibling_gate_pair_record_collection_sec",
        "sibling_gate_inflation_fit_sec",
        "sibling_gate_adjusted_tests_sec",
        "sibling_gate_fdr_sec",
        "traversal_sec",
    ):
        assert stage_timings[key] >= 0.0


def test_tbs_runner_accepts_explicit_fixed_sibling_gate_config() -> None:
    rng = np.random.default_rng(123)
    data = pd.DataFrame(
        rng.integers(0, 2, size=(16, 8)),
        index=[f"S{index}" for index in range(16)],
        columns=[f"F{index}" for index in range(8)],
    )
    result = run_tbs_on_distance(
        data,
        pdist(data.to_numpy(), metric="hamming"),
        DEFAULT_SIBLING_ALPHA,
        tree_linkage_method="average",
        sibling_gate_method="fixed_global_chi_square",
        sibling_gate_alpha_penalty=50.0,
        root_stability_guard_threshold=0.24,
        root_stability_subsample_replicates=12,
        root_stability_feature_fraction=0.8,
        root_stability_seed=0,
        root_selective_permutation_guard_replicates=1,
        root_selective_permutation_guard_seed=19,
        root_selective_permutation_guard_alpha=0.01,
        trace_level="full",
    )

    metadata = result.extra["gate_bundle"].metadata.config
    assert result.status == "ok"
    assert metadata.sibling_gate_profile_id is None
    assert metadata.sibling_gate_method == "fixed_global_chi_square"
    assert metadata.sibling_gate_alpha_penalty == 50.0
    assert metadata.root_stability_guard_threshold == 0.24
    assert metadata.root_stability_tree_distance_metric == "hamming"
    assert metadata.root_stability_tree_linkage_method == "average"
    assert metadata.root_selective_permutation_guard_replicates == 1
    assert metadata.root_selective_permutation_guard_seed == 19
    assert metadata.root_selective_permutation_guard_alpha == 0.01
    assert result.extra["sibling_gate_method"] == "fixed_global_chi_square"


def test_tbs_runner_selected_root_guard_blocks_known_categorical_false_root() -> None:
    from sklearn.metrics import adjusted_rand_score

    from benchmarks.diagnostics.calibration.sibling.gates.data_independent_sibling_gate_traversal_panel import (
        _generate_data_with_truth,
    )
    from benchmarks.validation.statistics.selected_edge_type1_geometry import case_contract
    from benchmarks.validation.sweeps.alpha_grid_search import select_cases

    case = next(case for case in select_cases(suite="full", case_names=("cat_clear_3cat_4c",)))
    (
        case_id,
        source_family,
        feature_representation,
        n_samples,
        n_features,
        n_categories,
    ) = case_contract(case)
    seed = 20309045
    results = {}
    for role in ("null", "signal"):
        data, feature_space, truth, _true_clusters = _generate_data_with_truth(
            case=case,
            case_id=case_id,
            source_family=source_family,
            feature_representation=feature_representation,
            n_samples=n_samples,
            n_features=n_features,
            n_categories=n_categories,
            data_role=role,
            seed=seed,
        )
        run = run_tbs_on_distance(
            data,
            pdist(data.to_numpy(dtype=float), metric="hamming"),
            0.01,
            tree_linkage_method="average",
            edge_alpha=0.001,
            feature_space=feature_space,
            sibling_gate_method="fixed_coordinate_bh",
            sibling_gate_alpha_penalty=50.0,
            root_stability_guard_threshold=0.24,
            root_stability_subsample_replicates=12,
            root_stability_feature_fraction=0.8,
            root_stability_seed=0,
            root_selective_permutation_guard_replicates=99,
            root_selective_permutation_guard_seed=0,
            root_selective_permutation_guard_alpha=0.01,
            root_selective_permutation_guard_scope="root",
            trace_level="full",
        )
        annotations = run.extra["annotations"]
        root = run.extra["tree"].root()
        metadata = run.extra["gate_bundle"].metadata.config
        assert metadata.sibling_gate_profile_id is None
        assert metadata.root_selective_permutation_guard_replicates == 99
        assert metadata.root_selective_permutation_guard_seed == 0
        assert metadata.root_selective_permutation_guard_alpha == 0.01
        results[role] = {
            "found_clusters": int(run.found_clusters),
            "ari": float(adjusted_rand_score(np.asarray(truth, dtype=int), np.asarray(run.labels))),
            "root_selective_p": float(annotations.at[root, "Selective_Permutation_P_Value"]),
            "blocked": bool(annotations.at[root, "Selective_Permutation_Guard_Blocked"]),
        }

    assert results["null"]["found_clusters"] == 1
    assert results["null"]["blocked"] is True
    assert results["null"]["root_selective_p"] == pytest.approx(0.17)
    assert results["signal"]["found_clusters"] == 4
    assert results["signal"]["blocked"] is False
    assert results["signal"]["root_selective_p"] == pytest.approx(0.01)
    assert results["signal"]["ari"] > 0.75
