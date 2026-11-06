"""A compact PC-stage surrogate for the cancer diffusion root failure."""

import numpy as np

from benchmarks.shared.generators.case_data_contracts import (
    CaseDataResult,
    case_metadata,
    continuous_dataframe_and_metadata,
    require_case_value,
)
from tree_break_selection.tree.continuous_distance import (
    CONTINUOUS_STANDARDIZED_EUCLIDEAN_TREE_DISTANCE_METRIC,
)


def generate_cancer_pc_surrogate(test_case: dict, seed: int | None) -> CaseDataResult:
    """Keep draws fixed across controls; change one row scale or two population means."""
    if tuple(test_case.get(k) for k in ("n_samples", "n_features", "n_clusters")) != (496, 30, 6):
        raise ValueError("Cancer PC surrogate requires n_samples=496, n_features=30, n_clusters=6.")
    hub_scale = float(require_case_value(test_case, "hub_scale", "Cancer PC surrogate"))
    isolation_offset = float(
        require_case_value(test_case, "isolation_offset", "Cancer PC surrogate")
    )
    for name, value in (("hub_scale", hub_scale), ("isolation_offset", isolation_offset)):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive.")

    sizes = [256, 64, 32, 16, 64, 64]
    labels = np.repeat(np.arange(6), sizes)
    means = np.zeros((6, 30))
    means[:4, 0] = np.arange(4) * 6.0
    means[4:, 1] = [-isolation_offset, isolation_offset]
    noise = np.random.default_rng(seed).normal(size=(496, 30))
    noise[0] *= hub_scale
    matrix = means[labels] + noise
    data, space, distances = continuous_dataframe_and_metadata(
        matrix,
        [f"S{i:04d}" for i in range(496)],
        [f"PC{i + 1:02d}" for i in range(30)],
        tree_distance_metric=CONTINUOUS_STANDARDIZED_EUCLIDEAN_TREE_DISTANCE_METRIC,
    )
    metadata = case_metadata(
        test_case=test_case,
        n_samples=496,
        n_features=30,
        n_clusters=6,
        noise=1.0,
        generator="cancer_pc_surrogate",
        source_family="cancer_pc_surrogate",
        feature_representation="continuous",
        requires_precomputed_tbs_distance=True,
        precomputed_distance_condensed=distances,
        distance_metric=CONTINUOUS_STANDARDIZED_EUCLIDEAN_TREE_DISTANCE_METRIC,
        extra={
            "feature_space": space,
            "simulation_model": "gaussian_mixture_with_one_low_variance_row",
            "benchmark_intent": "cancer_inspired_diffusion_root_failure_and_matched_controls",
            "scientific_caution": (
                "PC-stage mechanism surrogate; no gene counts or fitted cancer covariance. "
                "Rare populations are enlarged to preserve native k=15 connectivity. "
                "The hub has a different row covariance when hub_scale != 1. "
                "Population labels do not imply that every selected mean contrast is nonzero."
            ),
            "cluster_sizes": sizes,
            "population_means": means.tolist(),
            "hub_sample": "S0000",
            "hub_scale": hub_scale,
            "isolation_offset": isolation_offset,
            "noise_covariance": "identity_except_hub_scale_squared_times_identity_for_S0000",
            "source_diagnostic": "raw/inbox/choir-cell-lines/tbs-root-groups.md",
        },
    )
    return data, labels, matrix, metadata
