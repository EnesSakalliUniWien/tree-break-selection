from __future__ import annotations

import numpy as np
import pandas as pd


def selected_hierarchy_geometry_records(
    *,
    case_id: str = "case",
    case_category: str = "synthetic_family_a",
    source_family: str = "source_family_a",
    replicate_offset: int = 0,
) -> pd.DataFrame:
    """Build the shared synthetic geometry panel used by hierarchy studies."""
    rows: list[dict[str, object]] = []
    for index in range(20):
        replicate_index = replicate_offset + index // 2
        ratio = 1.0 + 0.2 * index
        cos2 = 0.05 + 0.04 * index
        rows.append(
            {
                "case_id": case_id,
                "case_category": case_category,
                "source_family": source_family,
                "feature_representation": "bernoulli",
                "feature_family": "bernoulli",
                "replicate_index": replicate_index,
                "selected_hierarchy_simulation_id": f"{case_id}:{replicate_index}",
                "feature_dimension": 40,
                "parent_sample_size": 20,
                "left_child_sample_size": 10,
                "right_child_sample_size": 10,
                "selected_hierarchy_ratio": ratio,
                "log_selected_hierarchy_ratio": float(np.log(ratio)),
                "parent_depth": index % 4,
                "parent_size_fraction": 1.0 - index / 40.0,
                "parent_size_bin": "root_0.75_1",
                "child_balance": 0.5 - index / 100.0,
                "child_size_ratio": 1.0 + index / 10.0,
                "branch_length_sum": 1.0 + index / 20.0,
                "branch_length_asymmetry": index / 50.0,
                "left_edge_raw_p_value": 0.001 + index / 1000.0,
                "right_edge_raw_p_value": 0.002 + index / 1000.0,
                "min_child_edge_raw_p_value": 0.001 + index / 1000.0,
                "max_child_edge_raw_p_value": 0.002 + index / 1000.0,
                "left_edge_bh_p_value": 0.01 + index / 1000.0,
                "right_edge_bh_p_value": 0.02 + index / 1000.0,
                "min_child_edge_bh_p_value": 0.01 + index / 1000.0,
                "max_child_edge_bh_p_value": 0.02 + index / 1000.0,
                "negative_log10_min_child_edge_bh_p_value": float(
                    -np.log10(0.01 + index / 1000.0)
                ),
                "raw_mp_signal_count": 2 + index % 3,
                "parent_test_projection_dimension": 2 + index % 3,
                "sibling_projection_dimension": 2 + index % 3,
                "effective_independent_rows": 20 + index,
                "mp_threshold_rows": 20 + index,
                "eigenvalue_effective_rank": 1.5 + index / 30.0,
                "top_eigenvalue_share": 0.8 - index / 100.0,
                "selected_eigenvalue_mass_fraction": 0.5 + index / 80.0,
                "eigengap_at_sibling_projection_dimension": 1.1 + index / 20.0,
                "selected_eigenvalue_over_mp_upper_bound": 0.7 + index / 50.0,
                "selected_subspace_cos2": cos2,
                "selected_subspace_sin2": 1.0 - cos2,
                "selected_subspace_tan2": (1.0 - cos2) / cos2,
                "top_component_cos2": 0.01 + (index % 5) / 30.0,
                "max_component_cos2": 0.02 + (index % 7) / 25.0,
                "component_cos2_entropy": 0.2 + ((index * 7) % 11) / 40.0,
            }
        )
    return pd.DataFrame.from_records(rows)
