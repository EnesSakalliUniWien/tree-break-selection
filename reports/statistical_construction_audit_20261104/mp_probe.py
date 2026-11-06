"""Isolate rank-selection and whitening behavior on fixed Gaussian leaf sets."""

import json
from pathlib import Path

import numpy as np

from tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald.projected_wald_kernel import (
    run_projected_wald_kernel,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.spectral.marchenko_pastur import (
    _process_node,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.spectral.node_spectral_task import (
    NodeSpectralTask,
)
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns


def solve(x, covariance, *, internal=()):
    m, d = x.shape
    space = continuous_feature_space_from_columns(tuple(range(d)))
    task = NodeSpectralTask(
        "node",
        tuple(range(m)),
        internal,
        x.mean(axis=0),
        space,
        {space.continuous_blocks[0].name: covariance},
    )
    return _process_node(task, x, 2, d, d)


def describe(result):
    return {
        "raw_mp_count": result.raw_mp_signal_count,
        "nominal_dimension": result.test_projection_dimension,
        "basis_rank": len(result.eigenvalues),
        "eigenvalue_min_max": (
            [float(result.eigenvalues[-1]), float(result.eigenvalues[0])]
            if len(result.eigenvalues)
            else []
        ),
    }


def main():
    rng = np.random.default_rng(20261104)
    m, d = 60, 30
    x = rng.standard_normal((m, d))
    covariance = np.cov(x, rowvar=False)
    known = solve(x, np.eye(d))
    fitted_full = solve(x, covariance)
    fitted_diagonal = solve(x, covariance.diagonal())
    x_small = rng.standard_normal((2, d))
    pair = solve(x_small, np.var(x_small, axis=0, ddof=1))
    zero = solve(np.zeros((2, d)), np.ones(d), internal=(np.zeros(d),))

    assert fitted_full.raw_mp_signal_count == 0
    assert fitted_full.test_projection_dimension == d
    assert np.max(np.abs(fitted_full.eigenvalues - (m - 1) / m)) < 1e-9
    assert pair.test_projection_dimension == 1
    assert np.isclose(pair.eigenvalues[0], d / 2)
    assert zero.mp_threshold_rows == 3 and zero.effective_independent_rows == 2

    def ranks_for(z):
        return [
            run_projected_wald_kernel(
                z,
                spectral_k=nominal,
                pca_projection=known.projection_matrix,
                pca_eigenvalues=known.eigenvalues,
                adaptive_dimension_energy_fraction=0.9,
            ).projection_dimension
            for nominal in (1, 2, 30)
        ]

    first = ranks_for(known.projection_matrix[0])
    last = ranks_for(known.projection_matrix[-1])
    assert first == [1, 1, 1] and last == [30, 30, 30]
    spike = np.zeros((m, d))
    spike[0, 0], spike[1, 0] = 20, -20
    assert solve(spike, np.eye(d)).raw_mp_signal_count == 1

    null_replays = 3000
    null_counts = {}
    for rows, cols in ((6, 30), (30, 30), (60, 30)):
        sample = rng.standard_normal((null_replays, rows, cols))
        sample -= sample.mean(axis=1, keepdims=True)
        sample = sample if rows < cols else sample.swapaxes(-1, -2)
        maximum = np.linalg.eigvalsh(sample @ sample.swapaxes(-1, -2) / rows)[:, -1]
        edges = {
            "implemented": (1 + np.sqrt(cols / rows)) ** 2,
            "centered_scale_matched": (np.sqrt(rows - 1) + np.sqrt(cols)) ** 2 / rows,
            "centered_unscaled": (1 + np.sqrt(cols / (rows - 1))) ** 2,
        }
        null_counts[f"{rows}x{cols}"] = {
            label: {"edge": float(edge), "count_above": int(np.sum(maximum > edge))}
            for label, edge in edges.items()
        }
    result = {
        "seed": 20261104,
        "null_replays_per_shape": null_replays,
        "known": describe(known),
        "fitted_full": describe(fitted_full),
        "fitted_diagonal": describe(fitted_diagonal),
        "fitted_diagonal_pair": describe(pair),
        "augmented_zero_rank": {
            "reported_threshold_rows": zero.mp_threshold_rows,
            "leaf_rows": zero.effective_independent_rows,
        },
        "adaptive_nominal_ranks": [1, 2, 30],
        "adaptive_aligned_first_ranks": first,
        "adaptive_aligned_last_ranks": last,
        "null_counts": null_counts,
    }
    serialized = json.dumps(result, indent=2) + "\n"
    destination = Path(__file__).parent / "generated/mp/probe.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(serialized)
    print(serialized, end="")


if __name__ == "__main__":
    main()
