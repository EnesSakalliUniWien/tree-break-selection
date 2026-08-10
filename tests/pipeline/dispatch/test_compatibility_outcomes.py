from __future__ import annotations

import numpy as np
import pandas as pd
from benchmarks.shared.runners.dispatch import run_clustering_result
from benchmarks.shared.runners.method_registry import METHOD_SPECS
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns

from .helpers import _toy_dataframe


def test_hamming_diffusion_dispatch_skips_continuous_feature_space() -> None:
    df = _toy_dataframe()

    result = run_clustering_result(
        data_df=df,
        method_id="tbs_diffusion",
        params=METHOD_SPECS["tbs_diffusion"].param_grid[0],
        seed=42,
        feature_space=continuous_feature_space_from_columns(df.columns),
    )

    assert result.status == "skip"
    assert result.labels is None
    assert result.found_clusters == 0
    assert "requires binary or one-hot" in str(result.skip_reason)
    assert result.extra == {
        "method_compatibility_status": "unsupported_feature_space",
        "method_compatibility_reason": "hamming_diffusion_requires_discrete_features",
    }

def test_adaptive_hamming_diffusion_dispatch_skips_continuous_feature_space() -> None:
    df = _toy_dataframe()

    result = run_clustering_result(
        data_df=df,
        method_id="tbs_diffusion_adaptive_nnls",
        params=METHOD_SPECS["tbs_diffusion_adaptive_nnls"].param_grid[0],
        seed=42,
        feature_space=continuous_feature_space_from_columns(df.columns),
    )

    assert result.status == "skip"
    assert result.labels is None
    assert "requires binary or one-hot" in str(result.skip_reason)


def test_adaptive_pydiffmap_dispatch_skips_zero_bandwidth_duplicate_blocks() -> None:
    values = np.vstack([np.zeros((8, 5)), np.eye(5)])
    df = pd.DataFrame(values, columns=[f"f{i}" for i in range(values.shape[1])])

    result = run_clustering_result(
        data_df=df,
        method_id="tbs_diffusion_adaptive_nnls",
        params=METHOD_SPECS["tbs_diffusion_adaptive_nnls"].param_grid[0],
        seed=42,
    )

    assert result.status == "skip"
    assert result.labels is None
    assert "zero local bandwidths" in str(result.skip_reason)
    assert result.extra == {
        "method_compatibility_status": "unsupported_duplicate_geometry",
        "method_compatibility_reason": "pydiffmap_variable_bandwidth_zero_bandwidth",
        "max_duplicate_count": 8,
        "pydiffmap_nnkde_query_size": 8,
    }
