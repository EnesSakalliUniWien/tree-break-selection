from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tree_break_selection.space_separation import (
    FEATURE_FAMILY_DIFFUSION_SETTINGS,
    feature_family_diffusion_geometry,
    feature_family_diffusion_unsupported_reason,
)
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns


def _frame(values: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame(
        values,
        index=[f"S{i}" for i in range(values.shape[0])],
        columns=[f"F{j}" for j in range(values.shape[1])],
    )


def _binary(n_samples: int = 60, n_features: int = 20, seed: int = 1) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    return _frame((rng.random((n_samples, n_features)) < 0.5).astype(float))


def test_continuous_data_gets_euclidean_diffusion_with_one_component_per_feature() -> None:
    values = np.random.default_rng(0).normal(size=(60, 8))
    values[30:] += 3.0
    data = _frame(values)
    feature_space = continuous_feature_space_from_columns(tuple(data.columns))

    geometry = feature_family_diffusion_geometry(data, feature_space)

    assert geometry.metadata["feature_family"] == "continuous"
    assert geometry.metadata["metric"] == "euclidean"
    assert geometry.metadata["diffusion_components"] == 8
    assert geometry.distance_condensed.shape == (60 * 59 // 2,)
    assert FEATURE_FAMILY_DIFFUSION_SETTINGS["continuous"].edge_projection_energy_fraction == 0.9


def test_binary_data_gets_hamming_diffusion_without_an_energy_fraction() -> None:
    geometry = feature_family_diffusion_geometry(_binary(), None)

    assert geometry.metadata["feature_family"] == "bernoulli"
    assert geometry.metadata["metric"] == "hamming"
    assert FEATURE_FAMILY_DIFFUSION_SETTINGS["bernoulli"].edge_projection_energy_fraction is None


def test_unsupported_reason_names_duplicate_rows_and_non_binary_values() -> None:
    clean = _binary()
    duplicated = clean.copy()
    duplicated.iloc[:8] = duplicated.iloc[0].to_numpy()

    assert feature_family_diffusion_unsupported_reason(clean, None) is None
    assert "identical rows" in feature_family_diffusion_unsupported_reason(duplicated, None)
    assert "binary" in feature_family_diffusion_unsupported_reason(clean * 2.0, None)
    with pytest.raises(ValueError, match="identical rows"):
        feature_family_diffusion_geometry(duplicated, None)
