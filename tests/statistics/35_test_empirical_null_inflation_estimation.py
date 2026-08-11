from __future__ import annotations

import math
from dataclasses import replace

import numpy as np
import pytest
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, f, gamma
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.inflation_correction.empirical_null_inflation_estimation import (
    CalibrationSupportThresholds,
    decide_empirical_null_calibration,
    decide_independent_unweighted_common_scale_calibration,
    fit_empirical_null_inflation_model,
    fit_independent_unweighted_common_scale_inflation_model,
    predict_empirical_inflation_factor,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.inflation_correction.inflation_adjusted_sibling_tests import (
    compute_inflation_adjusted_sibling_tests,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.inflation_correction.types.inflation_model import (
    CalibrationFitDiagnostics,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.types.sibling_pair_record import (
    SiblingPairRecord,
)


def _make_record(
    parent: str,
    *,
    stat: float,
    degrees_of_freedom: float,
    reference_scale: float = 1.0,
    is_null_like: bool = True,
    is_edge_blocked: bool = False,
    sibling_null_weight: float = 1.0,
    sibling_projection_dimension: float = 2.0,
    n_parent: int = 32,
    feature_family: str = "bernoulli",
    p_value: float = 0.5,
    calibration_dependency_group: object | None = None,
    calibration_dependency_weight: float | None = None,
) -> SiblingPairRecord:
    dependency_kwargs = {}
    if calibration_dependency_group is not None:
        dependency_kwargs["calibration_dependency_group"] = calibration_dependency_group
    if calibration_dependency_weight is not None:
        dependency_kwargs["calibration_dependency_weight"] = calibration_dependency_weight
    return SiblingPairRecord(
        parent=parent,
        left=f"{parent}L",
        right=f"{parent}R",
        stat=stat,
        reference_scale=reference_scale,
        degrees_of_freedom=degrees_of_freedom,
        p_value=p_value,
        branch_length_sum=0.1,
        n_parent=n_parent,
        is_null_like=is_null_like,
        is_edge_blocked=is_edge_blocked,
        sibling_null_weight=sibling_null_weight,
        sibling_projection_dimension=sibling_projection_dimension,
        feature_family=feature_family,
        **dependency_kwargs,
    )


@pytest.mark.parametrize(
    ("field_name", "invalid_value"),
    [
        ("min_supported_groups", 0),
        ("min_supported_groups", -1),
        ("min_supported_groups", True),
        ("min_supported_groups", 1.0),
        ("min_family_supported_groups", 0),
    ],
)
def test_calibration_support_thresholds_reject_invalid_counts(
    field_name: str,
    invalid_value: object,
) -> None:
    with pytest.raises(ValueError, match=field_name):
        CalibrationSupportThresholds(**{field_name: invalid_value})


@pytest.mark.parametrize(
    ("field_name", "invalid_value"),
    [
        ("min_family_effective_sample_size", 0.0),
        ("min_family_effective_sample_size", float("nan")),
        ("min_family_effective_sample_size", float("inf")),
        ("min_local_effective_sample_size", 0.999),
        ("min_local_effective_sample_size", float("-inf")),
    ],
)
def test_calibration_support_thresholds_reject_invalid_effective_sample_sizes(
    field_name: str,
    invalid_value: float,
) -> None:
    with pytest.raises(ValueError, match=field_name):
        CalibrationSupportThresholds(**{field_name: invalid_value})


@pytest.mark.parametrize(
    "invalid_value",
    [0.0, -0.1, 1.0001, float("nan"), float("inf")],
)
def test_calibration_support_thresholds_reject_invalid_max_weight_share(
    invalid_value: float,
) -> None:
    with pytest.raises(ValueError, match="max_weight_share"):
        CalibrationSupportThresholds(max_weight_share=invalid_value)


@pytest.mark.parametrize(
    "invalid_value",
    [-1e-12, float("nan"), float("inf")],
)
def test_calibration_support_thresholds_reject_invalid_leave_one_limit(
    invalid_value: float,
) -> None:
    with pytest.raises(ValueError, match="max_leave_one_group_delta_log_c"):
        CalibrationSupportThresholds(max_leave_one_group_delta_log_c=invalid_value)


def test_calibration_support_thresholds_accept_valid_boundaries() -> None:
    thresholds = CalibrationSupportThresholds(
        min_supported_groups=1,
        min_family_supported_groups=1,
        min_family_effective_sample_size=1.0,
        min_local_effective_sample_size=1.0,
        max_weight_share=1.0,
        max_leave_one_group_delta_log_c=0.0,
    )

    assert thresholds.max_leave_one_group_delta_log_c == 0.0


def test_fit_empirical_null_inflation_model_rejects_empty_calibration_set() -> None:
    with pytest.raises(ValueError, match="no sibling calibration records"):
        fit_empirical_null_inflation_model([])


def test_fit_empirical_null_inflation_model_rejects_nonfinite_statistic() -> None:
    records = [_make_record("bad_stat", stat=float("nan"), degrees_of_freedom=2.0)]

    with pytest.raises(ValueError, match="finite statistics"):
        fit_empirical_null_inflation_model(records)


def test_fit_empirical_null_inflation_model_rejects_no_positive_degrees_of_freedom() -> None:
    records = [_make_record("bad_df", stat=0.0, degrees_of_freedom=0.0)]

    with pytest.raises(ValueError, match="no positive-degree calibration records"):
        fit_empirical_null_inflation_model(records)


def test_fit_empirical_null_inflation_model_rejects_negative_statistic() -> None:
    records = [_make_record("negative", stat=-1.0, degrees_of_freedom=2.0)]

    with pytest.raises(ValueError, match="non-negative statistics"):
        fit_empirical_null_inflation_model(records)


def test_fit_empirical_null_inflation_model_rejects_invalid_null_weight() -> None:
    records = [
        _make_record(
            "bad_weight",
            stat=1.0,
            degrees_of_freedom=1.0,
            sibling_null_weight=1.5,
        )
    ]

    with pytest.raises(ValueError, match="sibling_null_weight"):
        fit_empirical_null_inflation_model(records)


def test_fit_empirical_null_inflation_model_excludes_selected_nonnull_records() -> None:
    records = [
        _make_record(
            "blocked",
            stat=4.0,
            degrees_of_freedom=1.0,
            is_edge_blocked=True,
            sibling_null_weight=0.5,
        ),
        _make_record(
            "focal",
            stat=20.0,
            degrees_of_freedom=4.0,
            is_null_like=False,
            sibling_null_weight=0.25,
        ),
    ]

    model = fit_empirical_null_inflation_model(records)

    assert model.n_calibration == 1
    assert model.n_edge_blocked_calibration == 1
    assert math.isclose(
        model.baseline_empirical_inflation_factor,
        4.0 / 1.0,
    )


def test_fit_empirical_null_inflation_model_reports_exclusive_support_roles() -> None:
    records = [
        _make_record(
            "tested_null",
            stat=2.0,
            degrees_of_freedom=1.0,
            is_null_like=True,
            is_edge_blocked=False,
        ),
        _make_record(
            "stopped",
            stat=8.0,
            degrees_of_freedom=2.0,
            is_null_like=True,
            is_edge_blocked=True,
        ),
    ]

    model = fit_empirical_null_inflation_model(records)

    assert model.n_calibration == 2
    assert model.n_strict_null_calibration == 1
    assert model.n_edge_blocked_calibration == 1
    assert model.sample_is_strict_null.tolist() == [True, False]
    assert model.sample_is_edge_blocked.tolist() == [False, True]


def test_fit_empirical_null_inflation_model_rejects_selected_nonnull_only_support() -> None:
    records = [
        _make_record(
            "selected",
            stat=20.0,
            degrees_of_freedom=4.0,
            is_null_like=False,
            sibling_null_weight=0.25,
        ),
    ]

    with pytest.raises(ValueError, match="selected non-null"):
        fit_empirical_null_inflation_model(records)


def _numeric_unit_scale_inflation_mle(records: list[SiblingPairRecord]) -> float:
    statistics = np.array([record.stat for record in records], dtype=float)
    degrees_of_freedom = np.array([record.degrees_of_freedom for record in records], dtype=float)
    weights = np.array([record.sibling_null_weight for record in records], dtype=float)

    def negative_log_likelihood(inflation_factor: float) -> float:
        if inflation_factor <= 0.0:
            return float("inf")
        return -float(
            np.sum(
                weights
                * gamma.logpdf(
                    statistics,
                    a=degrees_of_freedom / 2.0,
                    loc=0.0,
                    scale=2.0 * inflation_factor,
                )
            )
        )

    result = minimize_scalar(
        negative_log_likelihood,
        bounds=(1e-9, 100.0),
        method="bounded",
    )
    assert result.success
    return float(result.x)


def test_fit_empirical_null_inflation_model_matches_numeric_unit_scale_mle() -> None:
    records = [
        _make_record("p0", stat=2.0, degrees_of_freedom=1.0, sibling_null_weight=1.0),
        _make_record("p1", stat=9.0, degrees_of_freedom=3.0, sibling_null_weight=0.5),
        _make_record("p2", stat=20.0, degrees_of_freedom=5.0, sibling_null_weight=0.75),
    ]

    model = fit_empirical_null_inflation_model(records)
    numeric_mle = _numeric_unit_scale_inflation_mle(records)

    assert math.isclose(
        model.baseline_empirical_inflation_factor,
        numeric_mle,
        rel_tol=1e-6,
    )


def test_fit_empirical_null_inflation_model_rejects_non_unit_reference_scale() -> None:
    records = [
        _make_record("p0", stat=4.0, reference_scale=2.0, degrees_of_freedom=1.0),
        _make_record("p1", stat=12.0, reference_scale=3.0, degrees_of_freedom=2.0),
    ]

    with pytest.raises(ValueError, match="unit reference_scale"):
        fit_empirical_null_inflation_model(records)


def test_fit_empirical_null_inflation_model_uses_context_weighted_inflation() -> None:
    records = [
        _make_record("p0", stat=2.0, degrees_of_freedom=1.0, sibling_projection_dimension=1.0),
        _make_record("p1", stat=9.0, degrees_of_freedom=3.0, sibling_projection_dimension=2.0),
        _make_record("p2", stat=20.0, degrees_of_freedom=5.0, sibling_projection_dimension=4.0),
    ]

    model = fit_empirical_null_inflation_model(records)

    expected_scale = (2.0 + 9.0 + 20.0) / (1.0 + 3.0 + 5.0)

    assert math.isclose(
        model.baseline_empirical_inflation_factor,
        expected_scale,
        rel_tol=1e-9,
    )
    assert model.n_calibration == 3
    assert model.method == "context_weighted_supported_empirical_null_inflation"
    assert model.sample_statistics.tolist() == [2.0, 9.0, 20.0]
    assert model.sample_degrees_of_freedom.tolist() == [1.0, 3.0, 5.0]
    assert (
        predict_empirical_inflation_factor(
            model,
            records[0],
        )
        >= 1.0
    )


def test_selected_hierarchy_decision_preserves_diagnostics_but_fails_closed() -> None:
    records = [
        _make_record("strict", stat=4.0, degrees_of_freedom=2.0),
        _make_record(
            "edge_blocked",
            stat=8.0,
            degrees_of_freedom=2.0,
            is_null_like=True,
            is_edge_blocked=True,
        ),
    ]
    target = _make_record(
        "target",
        stat=16.0,
        degrees_of_freedom=2.0,
        is_null_like=False,
    )

    model = fit_empirical_null_inflation_model(records)
    decision = decide_empirical_null_calibration(model, target)

    assert decision.status == "undefined_unvalidated_reference_law"
    assert decision.c_hat is not None
    assert decision.p_value is None
    assert decision.support["n_supported_records"] == 2
    assert decision.support["n_positive_weight_records"] == 2
    assert decision.support["n_selected_nonnull_positive_weight_records"] == 0
    assert decision.support["n_strict_null_records"] == 1
    assert decision.support["n_edge_blocked_records"] == 1
    assert decision.support["n_family_supported_records"] == 2
    assert decision.support["support_contract_status"] == "below_internal_support_thresholds"
    assert "supported_groups_below_threshold" in str(
        decision.support["support_contract_failure_reasons"]
    )
    assert decision.exact_context["feature_family"] == "bernoulli"
    assert decision.estimator.endswith("family_baseline")
    assert decision.descriptive_strata["reason"] == (
        "same_selected_hierarchy_reference_law_unvalidated"
    )


def test_decide_empirical_null_calibration_can_enforce_support_thresholds() -> None:
    records = [
        _make_record("strict", stat=4.0, degrees_of_freedom=2.0),
        _make_record("blocked", stat=8.0, degrees_of_freedom=2.0, is_edge_blocked=True),
    ]
    target = _make_record(
        "target",
        stat=16.0,
        degrees_of_freedom=2.0,
        is_null_like=False,
    )

    model = fit_empirical_null_inflation_model(records)
    decision = decide_empirical_null_calibration(
        model,
        target,
        enforce_support_thresholds=True,
    )

    assert decision.status == "undefined_sparse_context"
    assert decision.c_hat is None
    assert decision.support["support_contract_status"] == "below_internal_support_thresholds"
    assert decision.descriptive_strata["reason"] == "internal_support_thresholds_failed"
    with pytest.raises(ValueError, match="undefined_sparse_context"):
        predict_empirical_inflation_factor(
            model,
            target,
            enforce_support_thresholds=True,
        )


def test_selected_hierarchy_decision_reports_passing_support_without_admissible_p_value() -> None:
    records = [
        _make_record("strict", stat=4.0, degrees_of_freedom=2.0),
        _make_record("blocked", stat=8.0, degrees_of_freedom=2.0, is_edge_blocked=True),
    ]
    target = _make_record(
        "target",
        stat=16.0,
        degrees_of_freedom=2.0,
        is_null_like=False,
    )
    permissive = CalibrationSupportThresholds(
        min_supported_groups=2,
        min_family_supported_groups=2,
        min_family_effective_sample_size=1.0,
        min_local_effective_sample_size=1.0,
        max_weight_share=1.0,
        max_leave_one_group_delta_log_c=10.0,
    )

    model = fit_empirical_null_inflation_model(records)
    decision = decide_empirical_null_calibration(
        model,
        target,
        enforce_support_thresholds=True,
        support_thresholds=permissive,
    )

    assert decision.status == "undefined_unvalidated_reference_law"
    assert decision.c_hat is not None
    assert decision.p_value is None
    assert decision.support["support_contract_status"] == ("passes_internal_support_thresholds")
    assert decision.support["support_contract_failure_reasons"] == ""


def test_compute_inflation_adjusted_sibling_tests_can_enforce_support_thresholds() -> None:
    records = [
        _make_record("strict", stat=4.0, degrees_of_freedom=2.0),
        _make_record("blocked", stat=8.0, degrees_of_freedom=2.0, is_edge_blocked=True),
        _make_record(
            "target",
            stat=16.0,
            degrees_of_freedom=2.0,
            is_null_like=False,
        ),
    ]
    model = fit_empirical_null_inflation_model(records)

    with pytest.raises(ValueError, match="undefined_unvalidated_reference_law"):
        compute_inflation_adjusted_sibling_tests(
            records,
            model=model,
        )

    with pytest.raises(ValueError, match="undefined_sparse_context"):
        compute_inflation_adjusted_sibling_tests(
            records,
            model=model,
            enforce_support_thresholds=True,
        )


def test_compute_inflation_adjusted_sibling_tests_rejects_missing_internal_model() -> None:
    records = [
        _make_record(
            "target",
            stat=160.0,
            degrees_of_freedom=2.0,
            is_null_like=False,
        ),
    ]

    with pytest.raises(ValueError, match="no internal model"):
        compute_inflation_adjusted_sibling_tests(records, model=None)


def test_decide_empirical_null_calibration_reports_missing_family_support() -> None:
    records = [
        _make_record("bernoulli", stat=4.0, degrees_of_freedom=2.0),
    ]
    target = _make_record(
        "target",
        stat=16.0,
        degrees_of_freedom=2.0,
        is_null_like=False,
        feature_family="categorical",
    )

    model = fit_empirical_null_inflation_model(records)
    decision = decide_empirical_null_calibration(model, target)

    assert decision.status == "undefined_no_family_support"
    assert decision.c_hat is None
    assert decision.p_value is None
    assert decision.support["n_supported_records"] == 1
    assert decision.support["n_family_supported_records"] == 0
    with pytest.raises(ValueError, match="undefined_no_family_support"):
        predict_empirical_inflation_factor(model, target)


def test_predict_empirical_inflation_conditions_on_parent_sample_size() -> None:
    records = [
        _make_record(
            "small_0",
            stat=80.0,
            degrees_of_freedom=2.0,
            sibling_projection_dimension=5.0,
            n_parent=2,
            feature_family="categorical",
        ),
        _make_record(
            "small_1",
            stat=60.0,
            degrees_of_freedom=2.0,
            sibling_projection_dimension=5.0,
            n_parent=3,
            feature_family="categorical",
        ),
        _make_record(
            "large_0",
            stat=10.0,
            degrees_of_freedom=5.0,
            sibling_projection_dimension=5.0,
            n_parent=100,
            feature_family="categorical",
        ),
    ]

    model = fit_empirical_null_inflation_model(records)

    small_parent_inflation = predict_empirical_inflation_factor(
        model,
        records[0],
    )
    large_parent_inflation = predict_empirical_inflation_factor(
        model,
        records[2],
    )

    assert large_parent_inflation < small_parent_inflation
    assert large_parent_inflation < model.baseline_empirical_inflation_factor


def test_fit_empirical_null_inflation_model_keeps_zero_ratios_as_calibration_data() -> None:
    records = [
        _make_record("zero", stat=0.0, degrees_of_freedom=2.0),
        _make_record("positive", stat=4.0, degrees_of_freedom=2.0),
    ]

    model = fit_empirical_null_inflation_model(records)

    assert model.n_calibration == 2
    assert model.baseline_empirical_inflation_factor == 1.0


def test_fit_empirical_null_inflation_model_enforces_one_sided_inflation_floor() -> None:
    records = [
        _make_record("p0", stat=0.4, degrees_of_freedom=1.0),
        _make_record("p1", stat=1.2, degrees_of_freedom=2.0),
    ]

    model = fit_empirical_null_inflation_model(records)

    assert model.baseline_empirical_inflation_factor == 1.0


def test_fit_empirical_null_inflation_model_reports_underflow_stable_effective_sample_size() -> (
    None
):
    records = [
        _make_record("p0", stat=4.0, degrees_of_freedom=1.0, sibling_null_weight=1e-240),
        _make_record("p1", stat=8.0, degrees_of_freedom=2.0, sibling_null_weight=2e-240),
    ]

    model = fit_empirical_null_inflation_model(records)

    assert math.isclose(model.effective_sample_size, 1.8)


@pytest.mark.parametrize(
    ("is_null_like", "is_edge_blocked", "expected"),
    [
        (True, False, True),
        (True, True, True),
        (False, False, False),
    ],
)
def test_sibling_pair_record_owns_empirical_null_support_rule(
    is_null_like: bool,
    is_edge_blocked: bool,
    expected: bool,
) -> None:
    record = _make_record(
        "p",
        stat=1.0,
        degrees_of_freedom=1.0,
        is_null_like=is_null_like,
        is_edge_blocked=is_edge_blocked,
    )

    assert record.has_empirical_null_support is expected


def test_sibling_pair_record_rejects_blocked_non_null_role() -> None:
    with pytest.raises(ValueError, match="is_edge_blocked implies is_null_like"):
        _make_record(
            "invalid_role",
            stat=1.0,
            degrees_of_freedom=1.0,
            is_null_like=False,
            is_edge_blocked=True,
        )


def test_fit_model_owns_one_aligned_immutable_calibration_sample() -> None:
    records = [
        _make_record("p0", stat=2.0, degrees_of_freedom=1.0),
        _make_record("p1", stat=8.0, degrees_of_freedom=2.0),
    ]

    sample = fit_empirical_null_inflation_model(records).sample

    assert sample.parent_ids == ("p0", "p1")
    assert sample.statistics.tolist() == [2.0, 8.0]
    assert sample.degrees_of_freedom.tolist() == [1.0, 2.0]
    assert sample.contexts.shape == (2, 2)
    assert sample.weights.shape == (2,)
    assert sample.statistics.flags.writeable is False
    with pytest.raises(ValueError, match="read-only"):
        sample.statistics[0] = 99.0


def test_fit_model_rejects_duplicate_calibration_parent_ids() -> None:
    records = [
        _make_record("duplicate", stat=2.0, degrees_of_freedom=1.0),
        _make_record("duplicate", stat=8.0, degrees_of_freedom=2.0),
    ]

    with pytest.raises(ValueError, match="parent_ids must be unique"):
        fit_empirical_null_inflation_model(records)


def test_fitted_model_rejects_inconsistent_or_nonfinite_summary_state() -> None:
    model = fit_empirical_null_inflation_model(
        [_make_record("cal", stat=2.0, degrees_of_freedom=2.0)]
    )

    with pytest.raises(ValueError, match="cannot be smaller than the calibration sample"):
        replace(
            model,
            fit_diagnostics=CalibrationFitDiagnostics(n_positive_weight_records=0),
        )
    with pytest.raises(ValueError, match="baseline_scale_estimate"):
        replace(model, baseline_scale_estimate=float("nan"))


def test_dependency_group_metrics_do_not_count_nested_stopping_records_as_independent() -> None:
    group_owner = _make_record(
        "stopping_group",
        stat=4.0,
        degrees_of_freedom=2.0,
        calibration_dependency_group="stopping_group",
        calibration_dependency_weight=1.0,
    )
    nested = _make_record(
        "nested_descendant",
        stat=100.0,
        degrees_of_freedom=2.0,
        is_edge_blocked=True,
        calibration_dependency_group="stopping_group",
        calibration_dependency_weight=1.0,
    )
    independent = _make_record(
        "independent_group",
        stat=4.0,
        degrees_of_freedom=2.0,
        calibration_dependency_group="independent_group",
        calibration_dependency_weight=1.0,
    )
    target = _make_record(
        "target",
        stat=8.0,
        degrees_of_freedom=2.0,
        is_null_like=False,
    )

    owner_only_decision = decide_empirical_null_calibration(
        fit_empirical_null_inflation_model([group_owner, independent]),
        target,
    )
    nested_decision = decide_empirical_null_calibration(
        fit_empirical_null_inflation_model([group_owner, nested, independent]),
        target,
    )

    assert owner_only_decision.support["n_supported_groups"] == 2
    assert nested_decision.support["n_supported_groups"] == 2
    assert owner_only_decision.support["family_effective_sample_size"] == pytest.approx(2.0)
    assert nested_decision.support["family_effective_sample_size"] == pytest.approx(2.0)
    assert owner_only_decision.support["local_max_group_weight_share"] == pytest.approx(0.5)
    assert nested_decision.support["local_max_group_weight_share"] == pytest.approx(0.5)
    assert nested_decision.support["leave_one_group_max_delta_log_c"] >= (
        owner_only_decision.support["leave_one_group_max_delta_log_c"]
    )


def test_independent_unweighted_common_scale_mode_uses_exact_f_reference_law() -> None:
    calibration_records = [
        _make_record(
            "cal0",
            stat=2.0,
            degrees_of_freedom=2.0,
            p_value=float(chi2.sf(2.0, 2.0)),
        ),
        _make_record(
            "cal1",
            stat=4.0,
            degrees_of_freedom=4.0,
            p_value=float(chi2.sf(4.0, 4.0)),
        ),
    ]
    target = _make_record(
        "target",
        stat=12.0,
        degrees_of_freedom=2.0,
        is_null_like=False,
        p_value=float(chi2.sf(12.0, 2.0)),
    )

    model = fit_independent_unweighted_common_scale_inflation_model(
        calibration_records,
        calibration_observation_ids_by_parent={
            "cal0": ("train0", "train1"),
            "cal1": ("train2",),
        },
    )
    decision = decide_independent_unweighted_common_scale_calibration(
        model,
        target,
        focal_observation_ids=("validation0", "validation1"),
    )

    expected = float(f.sf(6.0, 2.0, 6.0))
    assert decision.status == "independent_exact_admissible"
    assert decision.p_value == pytest.approx(expected)
    assert decision.c_hat == pytest.approx(1.0)
    assert decision.estimator == "exact_independent_unweighted_common_scale_f"


def test_exact_f_mode_accepts_one_shared_nonunit_reference_scale() -> None:
    calibration_records = [
        _make_record(
            "cal0",
            stat=4.0,
            reference_scale=2.0,
            degrees_of_freedom=2.0,
        ),
        _make_record(
            "cal1",
            stat=8.0,
            reference_scale=2.0,
            degrees_of_freedom=4.0,
        ),
    ]
    target = _make_record(
        "target",
        stat=24.0,
        reference_scale=2.0,
        degrees_of_freedom=2.0,
        is_null_like=False,
        p_value=float(chi2.sf(12.0, 2.0)),
    )

    model = fit_independent_unweighted_common_scale_inflation_model(
        calibration_records,
        calibration_observation_ids_by_parent={
            "cal0": ("train0",),
            "cal1": ("train1",),
        },
    )
    decision = decide_independent_unweighted_common_scale_calibration(
        model,
        target,
        focal_observation_ids=("validation",),
    )

    assert decision.p_value == pytest.approx(float(f.sf(6.0, 2.0, 6.0)))
    assert decision.exact_context["common_reference_scale"] == 2.0


def test_exact_f_mode_preserves_one_sided_unadjusted_p_value_floor() -> None:
    calibration_records = [
        _make_record(
            "cal0",
            stat=0.1,
            degrees_of_freedom=2.0,
            p_value=float(chi2.sf(0.1, 2.0)),
        ),
        _make_record(
            "cal1",
            stat=0.1,
            degrees_of_freedom=2.0,
            p_value=float(chi2.sf(0.1, 2.0)),
        ),
    ]
    raw_p_value = float(chi2.sf(4.0, 2.0))
    target = _make_record(
        "target",
        stat=4.0,
        degrees_of_freedom=2.0,
        is_null_like=False,
        p_value=raw_p_value,
    )

    model = fit_independent_unweighted_common_scale_inflation_model(
        calibration_records,
        calibration_observation_ids_by_parent={
            "cal0": ("train0",),
            "cal1": ("train1",),
        },
    )
    decision = decide_independent_unweighted_common_scale_calibration(
        model,
        target,
        focal_observation_ids=("validation",),
    )

    assert decision.c_hat == 1.0
    assert decision.p_value == pytest.approx(raw_p_value)


def test_exact_f_mode_rejects_weighted_or_overlapping_calibration_contracts() -> None:
    weighted = [
        _make_record(
            "cal",
            stat=2.0,
            degrees_of_freedom=2.0,
            sibling_null_weight=0.5,
        )
    ]
    with pytest.raises(ValueError, match="unit fixed weights"):
        fit_independent_unweighted_common_scale_inflation_model(
            weighted,
            calibration_observation_ids_by_parent={"cal": ("shared",)},
        )

    dependency_weighted = [
        _make_record(
            "cal",
            stat=2.0,
            degrees_of_freedom=2.0,
            calibration_dependency_group="cal",
            calibration_dependency_weight=0.5,
        )
    ]
    with pytest.raises(ValueError, match="unit fixed weights"):
        fit_independent_unweighted_common_scale_inflation_model(
            dependency_weighted,
            calibration_observation_ids_by_parent={"cal": ("shared",)},
        )

    calibration = [_make_record("cal", stat=2.0, degrees_of_freedom=2.0)]
    model = fit_independent_unweighted_common_scale_inflation_model(
        calibration,
        calibration_observation_ids_by_parent={"cal": ("shared",)},
    )
    target = _make_record(
        "target",
        stat=4.0,
        degrees_of_freedom=2.0,
        is_null_like=False,
        p_value=float(chi2.sf(4.0, 2.0)),
    )
    with pytest.raises(ValueError, match="disjoint"):
        decide_independent_unweighted_common_scale_calibration(
            model,
            target,
            focal_observation_ids=("shared",),
        )

    overlapping_calibration = [
        _make_record("cal0", stat=2.0, degrees_of_freedom=2.0),
        _make_record("cal1", stat=2.0, degrees_of_freedom=2.0),
    ]
    with pytest.raises(ValueError, match="pairwise disjoint"):
        fit_independent_unweighted_common_scale_inflation_model(
            overlapping_calibration,
            calibration_observation_ids_by_parent={
                "cal0": ("shared",),
                "cal1": ("shared",),
            },
        )

    unequal_scales = [
        _make_record("cal0", stat=2.0, reference_scale=1.0, degrees_of_freedom=2.0),
        _make_record("cal1", stat=4.0, reference_scale=2.0, degrees_of_freedom=2.0),
    ]
    with pytest.raises(ValueError, match="common reference_scale"):
        fit_independent_unweighted_common_scale_inflation_model(
            unequal_scales,
            calibration_observation_ids_by_parent={
                "cal0": ("train0",),
                "cal1": ("train1",),
            },
        )

    shared_dependency_group = [
        _make_record(
            "cal0",
            stat=2.0,
            degrees_of_freedom=2.0,
            calibration_dependency_group="shared_group",
            calibration_dependency_weight=1.0,
        ),
        _make_record(
            "cal1",
            stat=2.0,
            degrees_of_freedom=2.0,
            calibration_dependency_group="shared_group",
            calibration_dependency_weight=1.0,
        ),
    ]
    with pytest.raises(ValueError, match="distinct dependency groups"):
        fit_independent_unweighted_common_scale_inflation_model(
            shared_dependency_group,
            calibration_observation_ids_by_parent={
                "cal0": ("train0",),
                "cal1": ("train1",),
            },
        )


def test_zero_dimensional_selected_hierarchy_decision_remains_degenerate_not_unresolved() -> None:
    model = fit_empirical_null_inflation_model(
        [_make_record("cal", stat=2.0, degrees_of_freedom=2.0)]
    )
    target = _make_record(
        "zero",
        stat=0.0,
        degrees_of_freedom=0.0,
        is_null_like=False,
        sibling_projection_dimension=0.0,
        p_value=1.0,
    )

    decision = decide_empirical_null_calibration(model, target)

    assert decision.status == "internal_admissible"
    assert decision.p_value == 1.0
