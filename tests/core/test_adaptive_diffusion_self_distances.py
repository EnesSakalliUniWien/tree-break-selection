"""Identity masking before the discontinuous positive-distance epsilon statistic."""

import numpy as np
import pytest
from pydiffmap import kernel

from tree_break_selection.space_separation import diffusion
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns


@pytest.mark.parametrize(
    "entry_point",
    [
        "adaptive_diffusion_geometry",
        "feature_family_diffusion_geometry",
        "block_adaptive_diffusion_geometry",
    ],
)
def test_identity_mask_preserves_every_off_diagonal_value_and_explicit_zero(
    monkeypatch, entry_point
):
    original_fit = kernel.Kernel.fit
    original_resolver = diffusion.resolve_adaptive_epsilon
    original_coordinates = diffusion.compute_diffusion_coordinates
    before = []
    observed = []

    def fit(self, values):
        result = original_fit(self, values)
        # Deliberately emulate a positive self-distance artifact; do not rely on BLAS.
        self.scaled_dists.setdiag(1e-8)
        self.scaled_dists[0, 1] = 1e-12
        self.scaled_dists[1, 0] = 0.0
        before.append(self.scaled_dists.copy())
        return result

    def resolve(fitted, *args, **kwargs):
        graph = fitted.scaled_dists
        np.testing.assert_array_equal(graph.diagonal(), np.zeros(graph.shape[0]))
        difference = graph - before[0]
        difference.setdiag(0.0)
        assert np.count_nonzero(difference.data) == 0
        assert graph[0, 1] == 1e-12
        assert 0 in graph.getrow(1).indices  # Explicit off-diagonal zero survives.
        observed.append(True)
        return original_resolver(fitted, *args, **kwargs)

    def coordinates(weights, **kwargs):
        np.testing.assert_array_equal(weights.diagonal(), np.ones(weights.shape[0]))
        assert weights[1, 0] == 1.0  # The stored zero distance retains unit affinity.
        return original_coordinates(weights, **kwargs)

    monkeypatch.setattr(diffusion, "compute_diffusion_coordinates", coordinates)
    monkeypatch.setattr(kernel.Kernel, "fit", fit)
    monkeypatch.setattr(diffusion, "resolve_adaptive_epsilon", resolve)
    values = np.random.default_rng(4).normal(size=(20, 4))
    kwargs = {
        "k_neighbors": 15,
        "diffusion_time": 3,
        "n_components": 4,
        "metric": "euclidean",
        "bandwidth_type": "-1/(d+2)",
        "epsilon": "median",
    }
    if entry_point == "feature_family_diffusion_geometry":
        kwargs = {"feature_space": continuous_feature_space_from_columns(["a", "b", "c", "d"])}
    geometry = getattr(diffusion, entry_point)(values, **kwargs)
    assert observed == [True]
    assert (
        geometry.metadata["self_distance_policy"]
        == "zero_identity_diagonal_before_epsilon_and_kernel"
    )
