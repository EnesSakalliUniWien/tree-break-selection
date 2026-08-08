from __future__ import annotations

import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from benchmarks.diagnostics.calibration.selected.hierarchy.selected_hierarchy_geometry_covariates import (
    STUDY_ROLE,
    evaluate_candidate_equation_holdout,
    evaluate_candidate_equations,
    evaluate_covariate_block_models,
    evaluate_covariate_relationships,
    evaluate_selected_ratio_tail_law,
    run_selected_hierarchy_geometry_covariate_study,
    summarize_selected_geometry_by_case,
)
from tests.validation.calibration.selected.hierarchy.fixtures import (
    selected_hierarchy_geometry_records,
)


def _geometry_records() -> pd.DataFrame:
    return selected_hierarchy_geometry_records()


def _multi_case_geometry_records() -> pd.DataFrame:
    first = _geometry_records()
    second = _geometry_records().copy()
    second["case_id"] = "case_b"
    second["case_category"] = "synthetic_family_b"
    second["source_family"] = "source_family_b"
    second["replicate_index"] = second["replicate_index"] + 10
    second["selected_hierarchy_simulation_id"] = [
        f"case_b:{replicate_index}" for replicate_index in second["replicate_index"]
    ]
    second["selected_hierarchy_ratio"] = second["selected_hierarchy_ratio"] * 1.25
    second["log_selected_hierarchy_ratio"] = np.log(second["selected_hierarchy_ratio"])
    return pd.concat([first, second], ignore_index=True)


def test_relationships_report_descriptive_correlations() -> None:
    records = _geometry_records()

    relationships = evaluate_covariate_relationships(records)

    row = relationships[relationships["covariate"].eq("selected_subspace_cos2")].iloc[0]
    assert row["relationship_status"] == "descriptive_unadjusted"
    assert row["covariate_block"] == "angular"
    assert row["n_pairs"] == 20
    assert float(row["spearman_rho"]) > 0.9


def test_block_models_report_in_sample_descriptive_fit() -> None:
    records = _geometry_records()

    models = evaluate_covariate_block_models(records)

    angular = models[models["covariate_block"].eq("angular")].iloc[0]
    assert angular["model_status"] == "descriptive_in_sample"
    assert float(angular["r_squared"]) >= 0.0
    assert angular["study_role"] == STUDY_ROLE


def test_candidate_equations_score_mean_and_tail_behavior() -> None:
    records = _geometry_records()

    equations = evaluate_candidate_equations(records)

    edge = equations[equations["equation_id"].eq("edge_action")].iloc[0]
    assert edge["equation_status"] == "descriptive_in_sample"
    assert float(edge["mean_log_ratio_r_squared"]) >= 0.0
    assert 0.0 <= float(edge["tail_auc_from_linear_score"]) <= 1.0
    assert float(edge["tail_event_rate"]) > 0.0


def test_candidate_equations_report_replicate_and_case_holdout() -> None:
    records = _multi_case_geometry_records()

    holdout = evaluate_candidate_equation_holdout(records, n_replicate_folds=2)

    assert set(holdout["split_strategy"]) == {
        "replicate_modulo",
        "leave_one_case_out",
        "leave_one_source_family_out",
    }
    edge = holdout[
        holdout["equation_id"].eq("edge_action") & holdout["split_strategy"].eq("replicate_modulo")
    ].iloc[0]
    assert edge["equation_status"] == "descriptive_holdout_replicate_modulo"
    assert float(edge["n_test_rows"]) > 0
    assert 0.0 <= float(edge["holdout_tail_auc_from_linear_score"]) <= 1.0


def test_selected_ratio_tail_law_reports_context_support_and_holdout_error() -> None:
    records = _multi_case_geometry_records()

    tail_law = evaluate_selected_ratio_tail_law(
        records,
        n_folds=2,
        min_train_simulations=2,
        min_train_records=2,
        required_min_matching_simulations=100,
        required_min_matched_records=100,
    )

    assert "edge_action_bin" in tail_law.columns
    first = tail_law.iloc[0]
    assert first["tail_law_status"] == "descriptive_holdout_tail_law"
    assert not bool(first["production_tail_law_admissible"])
    assert "matching_simulations_below_tail_resolution_contract" in str(
        first["tail_law_admissibility_failure_reasons"]
    )
    assert float(first["heldout_exceedance_absolute_error"]) >= 0.0


def test_selected_ratio_tail_law_can_mark_supported_context_admissible() -> None:
    records = pd.concat([_geometry_records()] * 8, ignore_index=True)
    records["replicate_index"] = np.arange(records.shape[0])
    records["selected_hierarchy_simulation_id"] = [
        f"case:{replicate_index}" for replicate_index in records["replicate_index"]
    ]

    tail_law = evaluate_selected_ratio_tail_law(
        records,
        n_folds=4,
        min_train_simulations=10,
        min_train_records=10,
        required_min_matching_simulations=40,
        required_min_matched_records=40,
        max_exceedance_standard_error=0.2,
    )

    assert bool(tail_law.iloc[0]["production_tail_law_admissible"])
    assert tail_law.iloc[0]["tail_law_admissibility_failure_reasons"] == ""
    assert float(tail_law.iloc[0]["max_exceedance_standard_error"]) == 0.2


def test_selected_ratio_tail_law_counts_case_replicates_as_independent() -> None:
    first = _geometry_records()
    second = _geometry_records().copy()
    second["case_id"] = "case_b"
    second["source_family"] = first["source_family"].iloc[0]
    second["selected_hierarchy_simulation_id"] = [
        f"case_b:{replicate_index}" for replicate_index in second["replicate_index"]
    ]
    records = pd.concat([first, second], ignore_index=True)
    records["sibling_projection_dimension"] = 2
    records["negative_log10_min_child_edge_bh_p_value"] = 8.5

    tail_law = evaluate_selected_ratio_tail_law(
        records,
        n_folds=2,
        min_train_simulations=2,
        min_train_records=2,
        required_min_matching_simulations=15,
        required_min_matched_records=20,
        max_exceedance_standard_error=0.2,
    )

    assert tail_law.shape[0] == 1
    assert int(tail_law.iloc[0]["n_matching_simulations"]) == 20


def test_selected_ratio_tail_law_requires_independent_simulation_ids() -> None:
    records = _geometry_records().drop(columns=["selected_hierarchy_simulation_id"])

    with pytest.raises(KeyError, match="explicit independent simulation ids"):
        evaluate_selected_ratio_tail_law(records)


def test_case_summary_uses_selected_ratio_and_geometry_fields() -> None:
    summary = summarize_selected_geometry_by_case(_geometry_records())

    row = summary.iloc[0]
    assert row["case_id"] == "case"
    assert row["n_records"] == 20
    assert row["n_matching_simulations"] == 10
    assert row["selected_hierarchy_ratio_q95"] > row["selected_hierarchy_ratio_median"]
    assert row["selected_subspace_cos2_mean"] > 0.0


def test_smoke_run_writes_geometry_outputs() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        outputs = run_selected_hierarchy_geometry_covariate_study(
            case_names=["gauss_null_large"],
            output_dir=output_dir,
            n_replicates=2,
            seed=11,
        )

        assert "selected_geometry_records" not in outputs
        assert (output_dir / "case_summary.csv").exists()
        assert (output_dir / "geometry_summary_by_case.csv").exists()
        assert (output_dir / "covariate_relationships.csv").exists()
        assert (output_dir / "covariate_block_models.csv").exists()
        assert (output_dir / "candidate_equations.csv").exists()
        assert (output_dir / "candidate_equation_holdout.csv").exists()
        assert (output_dir / "selected_ratio_tail_law.csv").exists()
        assert (output_dir / "manifest.json").exists()


def test_smoke_run_can_write_optional_row_level_geometry() -> None:
    with tempfile.TemporaryDirectory() as tmp:
        output_dir = Path(tmp)
        outputs = run_selected_hierarchy_geometry_covariate_study(
            case_names=["gauss_null_large"],
            output_dir=output_dir,
            n_replicates=2,
            seed=11,
            write_selected_records=True,
        )

        assert "selected_geometry_records" in outputs
        assert (output_dir / "selected_geometry_records.csv").exists()


def test_relationships_require_declared_predictor_columns() -> None:
    with pytest.raises(KeyError, match="Missing predictor column"):
        evaluate_covariate_relationships(
            pd.DataFrame({"log_selected_hierarchy_ratio": [1.0, 2.0]}),
            covariates=(("angular", "selected_subspace_cos2"),),
        )
