"""Mixed discrete/continuous root selected-region law diagnostic.

The root selected-region replay separates smooth merge-margin evidence from
discrete tie-cell evidence. This panel joins the root margin summary and the
rank-aware tie-cell burden table into one diagnostic selected-region law table.

The result is diagnostic-only. It identifies which conditioning coordinates are
present and why root promotion must remain fail-closed until the discrete
tie-rank null law is calibrated.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from benchmarks.diagnostics.calibration.reporting import write_diagnostic_bundle
from benchmarks.diagnostics.calibration.root.root_tail_values import finite_float, require_columns

SCHEMA_VERSION = "root_selected_mixed_region_law/v2"
STUDY_ROLE = "diagnostic_root_selected_mixed_region_law_not_calibration"
GENERATED_BY = "benchmarks.diagnostics.calibration.root.selected.root_selected_mixed_region_law"

DEFAULT_MARGIN_TOLERANCE = 1e-12

ROWS_OUTPUT = "root_selected_mixed_region_law_rows.csv"
RELATIONSHIPS_OUTPUT = "root_selected_mixed_region_law_relationships.csv"
SUMMARY_OUTPUT = "root_selected_mixed_region_law_summary.csv"
MANIFEST_OUTPUT = "manifest.json"

ROW_COLUMNS = (
    "schema_version",
    "study_role",
    "case_id",
    "root_selected_region_law_status",
    "root_selected_region_event",
    "root_mixed_region_component",
    "root_calibration_status",
    "root_continuous_margin_status",
    "root_discrete_tie_status",
    "root_tie_rank_status",
    "root_sibling_selected_ratio",
    "root_child_balance",
    "root_child_construction_merge_count",
    "root_child_min_merge_margin",
    "root_child_min_merge_margin_abs",
    "root_child_near_active_merge_count",
    "root_child_tied_minimum_merge_count",
    "root_child_discrete_tie_cell_count",
    "root_child_smooth_constraint_count",
    "root_tie_step_count",
    "root_tie_step_fraction",
    "root_tie_cell_log_burden",
    "root_tie_cell_mean_log_multiplicity",
    "root_tie_cell_geometric_mean_multiplicity",
    "root_tie_rank_log_burden",
    "root_tie_rank_mean_fraction",
    "root_tie_rank_median_fraction",
    "root_tie_rank_to_tie_burden_fraction",
    "root_tie_rank_shortfall_log_burden",
    "root_edge_path_radial_distance",
    "root_edge_path_statistic_margin",
    "root_edge_extra_parent_projection_energy",
    "root_selected_eigenvalue_over_mp_upper_bound",
    "root_selected_eigenvalue_mass_fraction",
    "root_raw_mp_signal_count",
    "root_mp_threshold_rows",
    "root_active_feature_count",
    "root_full_eigenvalue_count",
    "root_full_component_eigenvalues_json",
    "root_projected_eigenvalues_json",
    "root_mp_upper_bound",
    "root_rank_fraction_edge_margin_product",
    "root_rank_fraction_spectral_product",
    "root_mixed_law_inference_status",
)

RELATIONSHIP_COLUMNS = (
    "schema_version",
    "study_role",
    "target",
    "covariate",
    "relationship_status",
    "valid_pair_count",
    "unique_covariate_count",
    "spearman_r",
    "spearman_p_value",
    "log_log_pearson_r",
)

SUMMARY_COLUMNS = (
    "schema_version",
    "study_role",
    "root_mixed_region_component",
    "root_calibration_status",
    "row_count",
    "median_root_sibling_selected_ratio",
    "median_tie_rank_fraction",
    "summary_status",
)

RELATIONSHIP_COVARIATES = (
    "root_child_balance",
    "root_tie_step_fraction",
    "root_tie_cell_log_burden",
    "root_tie_cell_mean_log_multiplicity",
    "root_tie_cell_geometric_mean_multiplicity",
    "root_tie_rank_log_burden",
    "root_tie_rank_mean_fraction",
    "root_tie_rank_median_fraction",
    "root_tie_rank_to_tie_burden_fraction",
    "root_tie_rank_shortfall_log_burden",
    "root_edge_path_radial_distance",
    "root_edge_path_statistic_margin",
    "root_edge_extra_parent_projection_energy",
    "root_selected_eigenvalue_over_mp_upper_bound",
    "root_selected_eigenvalue_mass_fraction",
    "root_raw_mp_signal_count",
    "root_rank_fraction_edge_margin_product",
    "root_rank_fraction_spectral_product",
)


@dataclass(frozen=True)
class RootSelectedMixedRegionLawConfig:
    """Input paths and numerical knobs for the mixed root-law panel."""

    output_dir: Path
    root_selected_region_summary_path: Path
    tie_cell_burden_rows_path: Path
    margin_tolerance: float = DEFAULT_MARGIN_TOLERANCE


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--root-selected-region-summary-path",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--tie-cell-burden-rows-path",
        type=Path,
        required=True,
    )
    parser.add_argument(
        "--margin-tolerance",
        type=float,
        default=DEFAULT_MARGIN_TOLERANCE,
    )
    return parser.parse_args()


def _finite_median(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce")
    numeric = numeric[np.isfinite(numeric)]
    return float(numeric.median()) if not numeric.empty else math.nan


def _safe_nonnegative_log1p(value: float) -> float:
    if not math.isfinite(value):
        return math.nan
    return float(np.log1p(max(value, 0.0)))


def _summary_lookup(frame: pd.DataFrame, label: str) -> dict[str, dict[str, object]]:
    if frame.empty:
        return {}
    require_columns(frame, {"case_id"}, label)
    return {str(row["case_id"]): dict(row) for _, row in frame.iterrows()}


def _region_component(
    *,
    law_status: str,
    discrete_count: int,
    tie_step_count: int,
    smooth_count: int,
) -> str:
    has_discrete = (
        law_status == "discrete_tie_cell_geometry_required"
        or discrete_count > 0
        or tie_step_count > 0
    )
    has_smooth = law_status == "smooth_first_order_signed_distance_defined" or smooth_count > 0
    if has_discrete and has_smooth:
        return "mixed_smooth_and_discrete_region"
    if has_discrete:
        return "discrete_tie_rank_region"
    if has_smooth:
        return "smooth_margin_region"
    return "undefined_root_selected_region"


def _continuous_margin_status(min_margin: float, tolerance: float) -> str:
    if not math.isfinite(min_margin):
        return "margin_missing"
    if abs(min_margin) <= tolerance:
        return "active_or_numerically_tied_margin"
    if min_margin > tolerance:
        return "positive_margin_observed"
    return "negative_margin_replay_warning"


def _discrete_tie_status(tie_step_count: int, tie_log_burden: float) -> str:
    if tie_step_count <= 0 and not (math.isfinite(tie_log_burden) and tie_log_burden > 0.0):
        return "no_discrete_tie_cell_observed"
    return "discrete_tie_cell_observed"


def _tie_rank_status(
    *,
    tie_step_count: int,
    rank_fraction: float,
    rank_log_burden: float,
) -> str:
    if tie_step_count <= 0:
        return "not_applicable_no_tie_cell"
    if math.isfinite(rank_fraction) and 0.0 < rank_fraction <= 1.0:
        return "selected_tie_rank_coordinate_observed"
    if math.isfinite(rank_log_burden):
        return "selected_tie_rank_log_observed_fraction_missing"
    return "selected_tie_rank_coordinate_missing"


def _calibration_status(component: str, tie_rank_status: str) -> str:
    if component in {
        "discrete_tie_rank_region",
        "mixed_smooth_and_discrete_region",
    }:
        if tie_rank_status.startswith("selected_tie_rank"):
            return "blocked_until_discrete_tie_rank_null_calibrated"
        return "blocked_until_discrete_tie_rank_coordinate_available"
    if component == "smooth_margin_region":
        return "smooth_first_order_candidate_diagnostic_only"
    return "unsupported_root_selected_region_geometry"


def _mixed_inference_status(
    *,
    component: str,
    calibration_status: str,
) -> str:
    if calibration_status.startswith("blocked_until_discrete"):
        return "fail_closed_discrete_root_law_missing"
    if component == "smooth_margin_region":
        return "smooth_margin_diagnostic_only_no_production_p_value"
    return "fail_closed_unsupported_root_region"


def build_root_selected_mixed_region_law_rows(
    *,
    root_summary: pd.DataFrame,
    tie_cell_burden_rows: pd.DataFrame,
    margin_tolerance: float = DEFAULT_MARGIN_TOLERANCE,
) -> pd.DataFrame:
    """Return one mixed root selected-region law row per case."""
    require_columns(
        root_summary,
        {
            "case_id",
            "root_selected_region_law_status",
            "root_sibling_selected_ratio",
            "root_child_balance",
        },
        "root selected-region summary",
    )
    require_columns(
        tie_cell_burden_rows,
        {
            "case_id",
            "root_tie_step_count",
            "root_tie_cell_log_burden",
            "root_tie_rank_log_burden",
            "root_tie_rank_median_fraction",
        },
        "tie-cell burden rows",
    )
    tie_by_case = _summary_lookup(tie_cell_burden_rows, "tie-cell burden rows")
    records: list[dict[str, object]] = []

    for _, root in root_summary.sort_values("case_id").iterrows():
        case_id = str(root["case_id"])
        tie = tie_by_case.get(case_id, {})
        law_status = str(root.get("root_selected_region_law_status", "missing"))
        min_margin = finite_float(root.get("root_child_min_merge_margin", math.nan))
        min_margin_abs = abs(min_margin) if math.isfinite(min_margin) else math.nan
        discrete_count = int(finite_float(root.get("root_child_discrete_tie_cell_count", 0)))
        smooth_count = int(finite_float(root.get("root_child_smooth_constraint_count", 0)))
        tie_step_count = int(finite_float(tie.get("root_tie_step_count", 0)))
        tie_log_burden = finite_float(tie.get("root_tie_cell_log_burden", math.nan))
        rank_log_burden = finite_float(tie.get("root_tie_rank_log_burden", math.nan))
        rank_mean_fraction = finite_float(tie.get("root_tie_rank_mean_fraction", math.nan))
        rank_median_fraction = finite_float(tie.get("root_tie_rank_median_fraction", math.nan))
        representative_rank_fraction = (
            rank_median_fraction if math.isfinite(rank_median_fraction) else rank_mean_fraction
        )
        rank_to_tie_fraction = (
            rank_log_burden / tie_log_burden
            if math.isfinite(rank_log_burden)
            and math.isfinite(tie_log_burden)
            and tie_log_burden > 0.0
            else math.nan
        )
        rank_shortfall = (
            max(tie_log_burden - rank_log_burden, 0.0)
            if math.isfinite(rank_log_burden) and math.isfinite(tie_log_burden)
            else math.nan
        )
        edge_margin = finite_float(root.get("root_edge_path_statistic_margin", math.nan))
        spectral_ratio = finite_float(
            root.get("root_selected_eigenvalue_over_mp_upper_bound", math.nan)
        )
        rank_edge_product = (
            representative_rank_fraction * _safe_nonnegative_log1p(edge_margin)
            if math.isfinite(representative_rank_fraction)
            else math.nan
        )
        rank_spectral_product = (
            representative_rank_fraction * max(spectral_ratio, 0.0)
            if math.isfinite(representative_rank_fraction) and math.isfinite(spectral_ratio)
            else math.nan
        )

        component = _region_component(
            law_status=law_status,
            discrete_count=discrete_count,
            tie_step_count=tie_step_count,
            smooth_count=smooth_count,
        )
        continuous_status = _continuous_margin_status(
            min_margin,
            float(margin_tolerance),
        )
        discrete_status = _discrete_tie_status(tie_step_count, tie_log_burden)
        rank_status = _tie_rank_status(
            tie_step_count=tie_step_count,
            rank_fraction=representative_rank_fraction,
            rank_log_burden=rank_log_burden,
        )
        calibration_status = _calibration_status(component, rank_status)
        inference_status = _mixed_inference_status(
            component=component,
            calibration_status=calibration_status,
        )

        records.append(
            {
                "schema_version": SCHEMA_VERSION,
                "study_role": STUDY_ROLE,
                "case_id": case_id,
                "root_selected_region_law_status": law_status,
                "root_selected_region_event": (
                    "E_root = E_margin intersect E_tie_cell intersect E_tie_rank"
                ),
                "root_mixed_region_component": component,
                "root_calibration_status": calibration_status,
                "root_continuous_margin_status": continuous_status,
                "root_discrete_tie_status": discrete_status,
                "root_tie_rank_status": rank_status,
                "root_sibling_selected_ratio": finite_float(
                    root.get("root_sibling_selected_ratio", math.nan)
                ),
                "root_child_balance": finite_float(root.get("root_child_balance", math.nan)),
                "root_child_construction_merge_count": int(
                    finite_float(root.get("root_child_construction_merge_count", 0))
                ),
                "root_child_min_merge_margin": min_margin,
                "root_child_min_merge_margin_abs": min_margin_abs,
                "root_child_near_active_merge_count": int(
                    finite_float(root.get("root_child_near_active_merge_count", 0))
                ),
                "root_child_tied_minimum_merge_count": int(
                    finite_float(root.get("root_child_tied_minimum_merge_count", 0))
                ),
                "root_child_discrete_tie_cell_count": discrete_count,
                "root_child_smooth_constraint_count": smooth_count,
                "root_tie_step_count": tie_step_count,
                "root_tie_step_fraction": finite_float(tie.get("root_tie_step_fraction", math.nan)),
                "root_tie_cell_log_burden": tie_log_burden,
                "root_tie_cell_mean_log_multiplicity": finite_float(
                    tie.get("root_tie_cell_mean_log_multiplicity", math.nan)
                ),
                "root_tie_cell_geometric_mean_multiplicity": finite_float(
                    tie.get(
                        "root_tie_cell_geometric_mean_multiplicity",
                        math.nan,
                    )
                ),
                "root_tie_rank_log_burden": rank_log_burden,
                "root_tie_rank_mean_fraction": rank_mean_fraction,
                "root_tie_rank_median_fraction": rank_median_fraction,
                "root_tie_rank_to_tie_burden_fraction": rank_to_tie_fraction,
                "root_tie_rank_shortfall_log_burden": rank_shortfall,
                "root_edge_path_radial_distance": finite_float(
                    root.get("root_edge_path_radial_distance", math.nan)
                ),
                "root_edge_path_statistic_margin": edge_margin,
                "root_edge_extra_parent_projection_energy": finite_float(
                    root.get("root_edge_extra_parent_projection_energy", math.nan)
                ),
                "root_selected_eigenvalue_over_mp_upper_bound": spectral_ratio,
                "root_selected_eigenvalue_mass_fraction": finite_float(
                    root.get("root_selected_eigenvalue_mass_fraction", math.nan)
                ),
                "root_raw_mp_signal_count": finite_float(
                    root.get("root_raw_mp_signal_count", math.nan)
                ),
                "root_mp_threshold_rows": finite_float(
                    root.get("root_mp_threshold_rows", math.nan)
                ),
                "root_active_feature_count": finite_float(
                    root.get("root_active_feature_count", math.nan)
                ),
                "root_full_eigenvalue_count": finite_float(
                    root.get("root_full_eigenvalue_count", math.nan)
                ),
                "root_full_component_eigenvalues_json": str(
                    root.get("root_full_component_eigenvalues_json", "")
                ),
                "root_projected_eigenvalues_json": str(
                    root.get("root_projected_eigenvalues_json", "")
                ),
                "root_mp_upper_bound": finite_float(root.get("root_mp_upper_bound", math.nan)),
                "root_rank_fraction_edge_margin_product": rank_edge_product,
                "root_rank_fraction_spectral_product": rank_spectral_product,
                "root_mixed_law_inference_status": inference_status,
            }
        )
    return pd.DataFrame.from_records(records, columns=ROW_COLUMNS)


def summarize_mixed_region_relationships(rows: pd.DataFrame) -> pd.DataFrame:
    """Relate mixed-law coordinates to root selected-ratio scale."""
    if rows.empty:
        return pd.DataFrame(columns=RELATIONSHIP_COLUMNS)
    require_columns(rows, {"root_sibling_selected_ratio"}, "mixed-law rows")
    target = pd.to_numeric(rows["root_sibling_selected_ratio"], errors="coerce")
    target_values = target.to_numpy(dtype=float)
    records: list[dict[str, object]] = []
    for covariate in RELATIONSHIP_COVARIATES:
        if covariate not in rows.columns:
            raise KeyError(f"Missing mixed-law covariate {covariate!r}.")
        values = pd.to_numeric(rows[covariate], errors="coerce")
        value_array = values.to_numpy(dtype=float)
        valid = np.isfinite(value_array) & np.isfinite(target_values)
        valid &= target_values > 0.0
        valid_count = int(np.count_nonzero(valid))
        unique_count = int(np.unique(value_array[valid]).shape[0])
        unique_target_count = int(np.unique(target_values[valid]).shape[0])
        if valid_count < 3:
            status = "insufficient_valid_pairs"
            spearman_r = math.nan
            spearman_p = math.nan
            log_log_pearson = math.nan
        elif unique_count < 2:
            status = "constant_covariate"
            spearman_r = math.nan
            spearman_p = math.nan
            log_log_pearson = math.nan
        elif unique_target_count < 2:
            status = "constant_target"
            spearman_r = math.nan
            spearman_p = math.nan
            log_log_pearson = math.nan
        else:
            status = "evaluated"
            spearman = stats.spearmanr(values[valid], target[valid])
            spearman_r = float(spearman.statistic)
            spearman_p = float(spearman.pvalue)
            positive_covariate = value_array[valid] > 0.0
            if bool(np.all(positive_covariate)):
                log_log_pearson = float(
                    np.corrcoef(
                        np.log(value_array[valid]),
                        np.log(target_values[valid]),
                    )[0, 1]
                )
            else:
                log_log_pearson = math.nan
        records.append(
            {
                "schema_version": SCHEMA_VERSION,
                "study_role": STUDY_ROLE,
                "target": "log_root_sibling_selected_ratio",
                "covariate": covariate,
                "relationship_status": status,
                "valid_pair_count": valid_count,
                "unique_covariate_count": unique_count,
                "spearman_r": spearman_r,
                "spearman_p_value": spearman_p,
                "log_log_pearson_r": log_log_pearson,
            }
        )
    return pd.DataFrame.from_records(records, columns=RELATIONSHIP_COLUMNS)


def summarize_mixed_region_law_rows(rows: pd.DataFrame) -> pd.DataFrame:
    """Summarize mixed root selected-region statuses."""
    if rows.empty:
        return pd.DataFrame(columns=SUMMARY_COLUMNS)
    records: list[dict[str, object]] = []
    group_columns = ["root_mixed_region_component", "root_calibration_status"]
    for (component, calibration_status), group in rows.groupby(
        group_columns,
        dropna=False,
        sort=True,
    ):
        records.append(
            {
                "schema_version": SCHEMA_VERSION,
                "study_role": STUDY_ROLE,
                "root_mixed_region_component": str(component),
                "root_calibration_status": str(calibration_status),
                "row_count": int(group.shape[0]),
                "median_root_sibling_selected_ratio": _finite_median(
                    group["root_sibling_selected_ratio"]
                ),
                "median_tie_rank_fraction": _finite_median(group["root_tie_rank_median_fraction"]),
                "summary_status": "diagnostic_only_not_calibrated",
            }
        )
    return pd.DataFrame.from_records(records, columns=SUMMARY_COLUMNS)


def evaluate_root_selected_mixed_region_law(
    config: RootSelectedMixedRegionLawConfig,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Read inputs and return mixed-law rows, relationships, and summary."""
    root_summary = pd.read_csv(config.root_selected_region_summary_path)
    tie_cell_burden_rows = pd.read_csv(config.tie_cell_burden_rows_path)
    rows = build_root_selected_mixed_region_law_rows(
        root_summary=root_summary,
        tie_cell_burden_rows=tie_cell_burden_rows,
        margin_tolerance=float(config.margin_tolerance),
    )
    relationships = summarize_mixed_region_relationships(rows)
    summary = summarize_mixed_region_law_rows(rows)
    return rows, relationships, summary


def run_root_selected_mixed_region_law(
    config: RootSelectedMixedRegionLawConfig,
) -> dict[str, Path]:
    """Write mixed root selected-region law outputs."""
    rows, relationships, summary = evaluate_root_selected_mixed_region_law(config)
    return write_diagnostic_bundle(
        output_dir=Path(config.output_dir),
        tables={"rows": rows, "relationships": relationships, "summary": summary},
        filenames={
            "rows": ROWS_OUTPUT,
            "relationships": RELATIONSHIPS_OUTPUT,
            "summary": SUMMARY_OUTPUT,
        },
        manifest={
            "schema_version": SCHEMA_VERSION,
            "study_role": STUDY_ROLE,
            "generated_by": GENERATED_BY,
            "config": config,
            "row_count": int(rows.shape[0]),
            "relationship_row_count": int(relationships.shape[0]),
            "summary_row_count": int(summary.shape[0]),
        },
        manifest_filename=MANIFEST_OUTPUT,
        include_row_counts=False,
    )


def main() -> None:
    args = parse_args()
    run_root_selected_mixed_region_law(
        RootSelectedMixedRegionLawConfig(
            output_dir=args.output_dir,
            root_selected_region_summary_path=args.root_selected_region_summary_path,
            tie_cell_burden_rows_path=args.tie_cell_burden_rows_path,
            margin_tolerance=float(args.margin_tolerance),
        )
    )


if __name__ == "__main__":
    main()
