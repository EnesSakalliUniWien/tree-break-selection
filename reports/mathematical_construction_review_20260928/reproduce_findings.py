"""Reproduce the 2026-09-28 review observations without changing method code.

Run from the repository root with its existing .venv/bin/python. Results describe
the source hashes recorded in findings.json; these are diagnostic probes, not a
full-pipeline calibration study or a new test suite.
"""

from __future__ import annotations

import hashlib
import json
import platform
import sys
from dataclasses import asdict
from importlib.metadata import version
from pathlib import Path

import matplotlib
import numpy as np
from scipy.linalg import subspace_angles
from scipy.spatial.distance import pdist
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from benchmarks.diagnostics.path_b.recursive_pvalue_geometry import principal_subspace_alignment
from tree_break_selection.hierarchy_analysis.decomposition.backends.eigen.decomposition import (
    eigendecompose_covariance,
)
from tree_break_selection.hierarchy_analysis.decomposition.backends.eigen.projection import (
    build_pca_projection,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald.projected_wald_kernel import (
    run_projected_wald_kernel,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.projection_dimension_estimation.projection_dimension_estimators import (
    estimate_marchenko_pastur_dimension,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.spectral.marchenko_pastur import (
    _process_node,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.spectral.node_spectral_task import (
    NodeSpectralTask,
)
from tree_break_selection.space_separation.diffusion import (
    compute_diffusion_coordinates,
    hamming_knn_diffusion_geometry,
)
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns


def adaptive_wald() -> dict:
    rotations = []
    for angle in (0.0, np.pi / 4):
        basis = np.array([[np.cos(angle), np.sin(angle)], [-np.sin(angle), np.cos(angle)]])
        result = run_projected_wald_kernel(
            np.array([3.0, 0.0]),
            spectral_k=2,
            pca_projection=basis,
            pca_eigenvalues=np.ones(2),
            adaptive_dimension_energy_fraction=0.9,
        )
        rotations.append({"angle_degrees": float(np.degrees(angle)), **asdict(result)})

    z = np.random.default_rng(1907).normal(size=(1_000_000, 2))
    energy = z * z
    total = energy.sum(axis=1)
    rank = 1 + (energy[:, 0] < 0.9 * total)
    statistic = np.where(rank == 1, energy[:, 0], total)
    adaptive_p = chi2.sf(statistic, rank)
    fixed_p = chi2.sf(total, 2)
    # Verify the vectorized experiment against the real kernel before using it.
    for index in range(128):
        actual = run_projected_wald_kernel(
            z[index],
            spectral_k=2,
            pca_projection=np.eye(2),
            pca_eigenvalues=np.ones(2),
            adaptive_dimension_energy_fraction=0.9,
        )
        assert actual.projection_dimension == rank[index]
        np.testing.assert_allclose(actual.p_value, adaptive_p[index], rtol=1e-12)
    rates = []
    for alpha in (0.001, 0.01, 0.05):
        count = int(np.count_nonzero(adaptive_p <= alpha))
        rate = count / len(z)
        rates.append(
            {
                "alpha": alpha,
                "adaptive_rejections": count,
                "adaptive_rate": rate,
                "monte_carlo_standard_error": float(np.sqrt(rate * (1.0 - rate) / len(z))),
                "fixed_rank_two_rejections": int(np.count_nonzero(fixed_p <= alpha)),
            }
        )
    return {
        "rotations": rotations,
        "seed": 1907,
        "independent_contrasts": len(z),
        "kernel_cross_checks": 128,
        "null_rates": rates,
    }


def diffusion_probes() -> dict:
    weights = np.zeros((6, 6))
    for index in range(6):
        weights[index, (index + 1) % 6] = weights[(index + 1) % 6, index] = 1.0
    permutation = np.array([2, 5, 0, 4, 1, 3])
    errors = {}
    for rank in (1, 2):
        before = compute_diffusion_coordinates(weights, diffusion_time=2, n_components=rank)
        after = compute_diffusion_coordinates(
            weights[np.ix_(permutation, permutation)],
            diffusion_time=2,
            n_components=rank,
        )[np.argsort(permutation)]
        errors[str(rank)] = float(np.max(np.abs(pdist(before) - pdist(after))))
    binary = np.array([[1, 0, 0], [1, 1, 0], [0, 1, 1], [0, 0, 1]], dtype=float)
    hamming = hamming_knn_diffusion_geometry(
        binary, k_neighbors=1, diffusion_time=2, n_components=3
    )
    path = np.array([[0.0, 1.0, 0.0], [1.0, 0.0, 1.0], [0.0, 1.0, 0.0]])
    markov = path / path.sum(axis=1)[:, None]
    path_coordinates = compute_diffusion_coordinates(path, diffusion_time=2, n_components=2)
    return {
        "six_cycle_normalized_spectrum": np.linalg.eigvalsh(weights / 2.0).tolist(),
        "permutation": permutation.tolist(),
        "permutation_distance_errors": errors,
        "hamming_k_one_distances": hamming.distance_condensed.tolist(),
        "path_three_clipped_symmetric_distances": pdist(path_coordinates).tolist(),
        "path_three_markov_power_row_distances": pdist(markov @ markov).tolist(),
    }


def pca_probes() -> dict:
    data = np.array([[1.0, 0.0], [-1.0, 0.0], [0.0, 1.0], [0.0, -1.0]])
    cases = []
    for permutation in (np.array([0, 1]), np.array([1, 0])):
        eig = eigendecompose_covariance(data[:, permutation], compute_eigenvectors=True)
        dimension = estimate_marchenko_pastur_dimension(eig.eigenvalues, n_samples=4, n_features=2)
        basis, _ = build_pca_projection(
            eig,
            projection_dimension=dimension.test_projection_dimension,
            n_features_total=2,
        )
        basis = basis[:, np.argsort(permutation)]
        cases.append(
            {
                "feature_permutation": permutation.tolist(),
                "eigenvalues": eig.eigenvalues.tolist(),
                "dimension": asdict(dimension),
                "physical_projector": (basis.T @ basis).tolist(),
            }
        )
    difference = np.linalg.norm(
        np.asarray(cases[0]["physical_projector"]) - cases[1]["physical_projector"],
    )
    base = np.random.default_rng(104).normal(size=(20, 2))
    space = continuous_feature_space_from_columns(["x", "y"])
    block_name = space.continuous_blocks[0].name
    whitening = []
    for separation in (0.0, 10.0):
        values = base.copy()
        values[:10, 0] -= separation
        values[10:, 0] += separation
        covariance = np.cov(values, rowvar=False, ddof=1)
        task = NodeSpectralTask(
            node_id="root",
            row_indices=tuple(range(20)),
            internal_distributions=(),
            null_distribution=values.mean(axis=0),
            feature_space=space,
            continuous_covariance_by_block={block_name: covariance},
        )
        result = _process_node(task, values, 1, None, 2, True)
        whitening.append(
            {
                "separation": separation,
                "raw_eigenvalues": np.linalg.eigvalsh(covariance).tolist(),
                "whitened_eigenvalues": result.full_eigenvalues.tolist(),
                "raw_mp_signal_count": result.raw_mp_signal_count,
                "test_projection_dimension": result.test_projection_dimension,
            }
        )
    return {
        "tied_covariance": cases,
        "projector_difference_frobenius": float(difference),
        "continuous_whitening_seed": 104,
        "continuous_whitening": whitening,
    }


def principal_angle_probe() -> dict:
    parent, child = np.eye(2), np.array([[0.0, 1.0]])
    return {
        "original_basis": principal_subspace_alignment(parent, child),
        "reordered_parent_basis": principal_subspace_alignment(parent[::-1], child),
        "full_rectangular_cross_gram_singular_values": np.linalg.svd(
            parent @ child.T, compute_uv=False
        ).tolist(),
        "scipy_full_subspace_angles_degrees": np.degrees(
            subspace_angles(parent.T, child.T)
        ).tolist(),
    }


def line_angles(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    left = left / np.linalg.norm(left, axis=0)
    right = right / np.linalg.norm(right, axis=0)
    return np.degrees(np.arccos(np.clip(np.abs(left.T @ right), 0.0, 1.0)))


def paired_eigenvectors() -> dict:
    weights = np.array([[0, 1, 0, 0], [1, 0, 2, 0], [0, 2, 0, 3], [0, 0, 3, 0]], dtype=float)
    degrees = weights.sum(axis=1)
    markov = weights / degrees[:, None]
    symmetric = weights / np.sqrt(degrees[:, None] * degrees[None, :])
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    right = eigenvectors / np.sqrt(degrees[:, None])
    left = eigenvectors * np.sqrt(degrees[:, None])
    panels = [line_angles(right, right), line_angles(left, left), line_angles(left, right)]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2), layout="constrained")
    titles = ["Right with right", "Left with left", "Left with right"]
    for ax, angles, title in zip(axes, panels, titles, strict=True):
        graphic = ax.imshow(angles, vmin=0, vmax=90, cmap="cividis")
        ax.set_title(title)
        ax.set_xticks(range(4), labels=["1", "2", "3", "4"])
        ax.set_yticks(range(4), labels=["1", "2", "3", "4"])
        ax.set_xlabel("Eigenvector index")
        for row in range(4):
            for column in range(4):
                ax.text(
                    column,
                    row,
                    f"{angles[row, column]:.1f}°",
                    ha="center",
                    va="center",
                    color="white" if angles[row, column] < 45 else "black",
                    fontsize=10,
                )
    axes[0].set_ylabel("Eigenvector index")
    fig.colorbar(graphic, ax=axes, label="Acute angle between vector lines (degrees)", shrink=0.8)
    fig.suptitle("Paired eigenvectors of P = D⁻¹W: Euclidean angles", fontsize=14)
    fig.savefig(OUTPUT / "paired_eigenvector_angles.png", dpi=170)
    fig.savefig(OUTPUT / "paired_eigenvector_angles.pdf")
    plt.close(fig)
    return {
        "illustrative_operator_not_benchmark_data": True,
        "weights": weights.tolist(),
        "degrees": degrees.tolist(),
        "eigenvalues_ascending": eigenvalues.tolist(),
        "right_vectors": right.tolist(),
        "left_vectors": left.tolist(),
        "left_right_gram": (left.T @ right).tolist(),
        "left_left_gram": (left.T @ left).tolist(),
        "right_right_gram": (right.T @ right).tolist(),
        "right_right_line_angles_degrees": panels[0].tolist(),
        "left_left_line_angles_degrees": panels[1].tolist(),
        "left_right_line_angles_degrees": panels[2].tolist(),
        "right_residual": float(np.linalg.norm(markov @ right - right * eigenvalues)),
        "left_residual": float(np.linalg.norm(markov.T @ left - left * eigenvalues)),
        "biorthogonality_residual": float(np.linalg.norm(left.T @ right - np.eye(4))),
        "right_degree_metric_residual": float(
            np.linalg.norm(right.T @ (degrees[:, None] * right) - np.eye(4))
        ),
        "left_inverse_degree_metric_residual": float(
            np.linalg.norm(left.T @ (left / degrees[:, None]) - np.eye(4))
        ),
        "paired_middle_subspaces_angles_degrees": np.degrees(
            subspace_angles(left[:, 1:3], right[:, 1:3])
        ).tolist(),
    }


def main() -> None:
    modules = [
        sys.modules[run_projected_wald_kernel.__module__],
        sys.modules[compute_diffusion_coordinates.__module__],
        sys.modules[eigendecompose_covariance.__module__],
        sys.modules[build_pca_projection.__module__],
        sys.modules[estimate_marchenko_pastur_dimension.__module__],
        sys.modules[_process_node.__module__],
        sys.modules[principal_subspace_alignment.__module__],
    ]
    source_paths = [Path(module.__file__).resolve() for module in modules]
    source_paths.extend(
        ROOT / name
        for name in [
            "README.md",
            "pyproject.toml",
            "uv.lock",
            "Makefile",
            "tree_break_selection/hierarchy_analysis/statistics/contrast_covariance.py",
            "tree_break_selection/tree/distributions.py",
            "tree_break_selection/space_separation/adaptive_cosine.py",
            "tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_reference_distribution.py",
            "tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_projection_basis.py",
        ]
    )
    source_paths.append(Path(__file__).resolve())
    result = {
        "review_date": "2026-09-28",
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {
            name: version(name)
            for name in (
                "numpy",
                "scipy",
                "pandas",
                "matplotlib",
                "pytest",
                "scikit-learn",
                "pydiffmap",
            )
        },
        "source_sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in source_paths
        },
        "adaptive_wald": adaptive_wald(),
        "diffusion": diffusion_probes(),
        "pca": pca_probes(),
        "principal_angles": principal_angle_probe(),
        "paired_eigenvectors": paired_eigenvectors(),
    }
    (OUTPUT / "findings.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            {
                "output": str(OUTPUT),
                "null_contrasts": 1_000_000,
                "kernel_cross_checks": 128,
                "source_hashes": len(result["source_sha256"]),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
