"""Discovery/test separation for paired independent Gaussian measurements.

Fit topology and PCA using discovery measurements only. Test a separate matrix
with the SAME leaf identities and a known common row covariance. Observation IDs
must distinguish measurements, not leaves. These checks enforce provenance and
reject obvious reuse; independence and the known Gaussian model are scientific
input contracts, not properties inferred from labels or observed values.

This explicit API cannot replace an ordinary single-matrix call: that call lacks
independent test measurements and known covariance. No empirical inflation is
fitted here. Existing Tree-BH, sibling BH and traversal are reused.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass

import networkx as nx
import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.spatial.distance import pdist

from tree_break_selection.eigenvalue_groups import tied_eigenvalue_neighbors
from tree_break_selection.hierarchy_analysis.decomposition.backends.eigen.decomposition import (
    eigendecompose_covariance,
)
from tree_break_selection.hierarchy_analysis.decomposition.backends.eigen.projection import (
    build_pca_projection,
)
from tree_break_selection.hierarchy_analysis.statistics.multiple_testing import (
    apply_tree_bh_correction,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald.energy_selected_reference import (
    energy_selected_gaussian_sf,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald.projected_wald_kernel import (
    run_projected_wald_kernel,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.inflated_projected_wald_annotation.fdr_annotation import (
    apply_traversal_aligned_sibling_bh_results,
    init_sibling_annotation_df,
)
from tree_break_selection.hierarchy_analysis.tree_decomposition import TreeDecomposition
from tree_break_selection.space_separation import adaptive_diffusion_geometry
from tree_break_selection.tree.construction import DEFAULT_TREE_LINKAGE_METHOD, build_tree
from tree_break_selection.tree.optimized_branch_lengths import (
    fit_fixed_topology_nnls_branch_lengths,
)
from tree_break_selection.tree.poset_tree import PosetTree


def _readonly(values):
    array = np.array(values, dtype=float, copy=True)
    array.setflags(write=False)
    return array


def _digest(values):
    return hashlib.sha256(np.ascontiguousarray(values, dtype=float).tobytes()).hexdigest()


def _observations(ids, n):
    ids = tuple(ids)
    if len(ids) != n or len(set(ids)) != n or any(x is None for x in ids):
        raise ValueError("Observation IDs must be non-null, unique and aligned with rows.")
    return frozenset(ids)


def _values(data):
    if not isinstance(data, pd.DataFrame) or not data.index.is_unique or not data.columns.is_unique:
        raise ValueError("Measurements require a DataFrame with unique leaf and feature labels.")
    if any(not isinstance(label, str) for label in data.index):
        raise ValueError("Leaf identities must be strings for all supported tree builders.")
    values = data.to_numpy(dtype=float, copy=True)
    if values.ndim != 2 or min(values.shape) < 1 or not np.isfinite(values).all():
        raise ValueError("Measurements must be nonempty and finite.")
    return values


@dataclass(frozen=True)
class IndependentGaussianNode:
    parent: object
    left: object
    right: object
    left_rows: tuple[int, ...]
    right_rows: tuple[int, ...]
    projection: np.ndarray
    eigenvalues: np.ndarray


@dataclass(frozen=True)
class IndependentGaussianDiffusion:
    """Discovery-only Euclidean pydiffmap and optional fixed-topology NNLS.

    Defaults match the adaptive NNLS runner except that continuous Gaussian
    discovery uses Euclidean distances instead of its binary Hamming preset.
    NNLS targets original discovery columns, not diffusion coordinates.
    Fitted lengths describe geometry, not the known test covariance.
    """

    k_neighbors: int = 10
    diffusion_time: int = 3
    n_components: int = 30
    bandwidth_type: str | float | None = "-1/(d+2)"
    epsilon: str | float = "median"
    fit_nnls: bool = True
    nnls_pair_sample_size: int | None = 50_000
    nnls_random_state: int = 0
    nnls_solver_tolerance: float = 1e-5
    nnls_max_iterations: int | None = 1000


@dataclass(frozen=True)
class IndependentGaussianTree:
    tree: PosetTree
    leaf_labels: tuple[str, ...]
    feature_labels: tuple[object, ...]
    discovery_observation_ids: frozenset[object]
    discovery_sha256: str
    covariance_cholesky: np.ndarray
    nodes: tuple[IndependentGaussianNode, ...]
    energy_fraction: float | None
    tree_edges: frozenset[tuple[object, object]]
    tree_branch_lengths: tuple[tuple[object, object, float], ...]
    diffusion_coordinates: np.ndarray | None
    geometry_metadata_json: str
    branch_length_metadata_json: str


@dataclass(frozen=True)
class IndependentGaussianResult:
    tests: pd.DataFrame
    annotations: pd.DataFrame
    decomposition: dict
    inference_sha256: str
    reference_assumptions: str = "independent_rows_known_common_gaussian_covariance"


def fit_independent_gaussian_tree(
    discovery: pd.DataFrame,
    *,
    known_covariance: np.ndarray,
    observation_ids,
    linkage_method: str = DEFAULT_TREE_LINKAGE_METHOD,
    energy_fraction: float | None = 0.9,
    diffusion: IndependentGaussianDiffusion | None = None,
) -> IndependentGaussianTree:
    """Freeze topology and positive PCA candidate directions on discovery data.

    Each discovery/test row is a separate Gaussian measurement of the same leaf
    mean, with the supplied common covariance. Rows and the two measurement
    matrices must be mutually independent. Covariance must be known externally,
    not estimated from either matrix and relabeled as known.
    """
    values = _values(discovery)
    ids = _observations(observation_ids, len(values))
    if len(values) < 2:
        raise ValueError("Independent tree fitting requires at least two leaves.")
    if energy_fraction is not None and not 0 < float(energy_fraction) <= 1:
        raise ValueError(
            "Energy fraction must be in (0, 1], or None for fixed full candidate rank."
        )
    covariance = np.asarray(known_covariance, dtype=float)
    if covariance.shape != (values.shape[1], values.shape[1]) or not np.isfinite(covariance).all():
        raise ValueError("Known covariance must be a finite square feature covariance.")
    if not np.allclose(covariance, covariance.T, atol=1e-12, rtol=1e-12):
        raise ValueError("Known covariance must be symmetric.")
    try:
        cholesky = np.linalg.cholesky(covariance)
    except np.linalg.LinAlgError as exc:
        raise ValueError("Known covariance must be positive definite.") from exc
    whitened = solve_triangular(cholesky, values.T, lower=True).T
    diffusion_coordinates = None
    geometry_metadata = {"method": "known_covariance_euclidean", "input": "discovery_only"}
    if diffusion is None:
        distances = pdist(whitened)
    else:
        geometry = adaptive_diffusion_geometry(
            discovery,
            k_neighbors=diffusion.k_neighbors,
            diffusion_time=diffusion.diffusion_time,
            n_components=diffusion.n_components,
            metric="euclidean",
            bandwidth_type=diffusion.bandwidth_type,
            epsilon=diffusion.epsilon,
        )
        distances = geometry.distance_condensed
        diffusion_coordinates = _readonly(geometry.coordinates)
        geometry_metadata = {
            **geometry.metadata,
            "method": "adaptive_pydiffmap_diffusion",
            "input": "original_discovery_measurements_only",
            "coordinates_sha256": _digest(geometry.coordinates),
            "distances_sha256": _digest(distances),
        }
    tree = build_tree(
        discovery,
        distances,
        linkage_method=linkage_method,
    ).tree
    branch_metadata = {"method": "constructor_lengths", "covariance_role": "none"}
    if diffusion is not None and diffusion.fit_nnls:
        optimization = fit_fixed_topology_nnls_branch_lengths(
            tree,
            discovery,
            target_metric="mean_squared_standardized_euclidean",
            pair_sample_size=diffusion.nnls_pair_sample_size,
            random_state=diffusion.nnls_random_state,
            solver_tolerance=diffusion.nnls_solver_tolerance,
            max_iterations=diffusion.nnls_max_iterations,
            apply_nonconverged=False,
        )
        if optimization.status != "ok" or not optimization.applied_to_tree:
            raise ValueError(
                "Discovery NNLS did not converge: " + json.dumps(optimization.to_dict())
            )
        branch_metadata = {
            **optimization.to_dict(),
            "input": "original_discovery_measurements_only",
            "covariance_role": "none",
        }
    descendants = tree.compute_descendant_sets(use_labels=True)
    positions = {name: i for i, name in enumerate(discovery.index)}
    nodes = []
    for parent in tree.nodes:
        children = tuple(tree.successors(parent))
        if not children:
            continue
        if len(children) != 2:
            raise ValueError("Independent Gaussian gates require a binary discovery tree.")
        left, right = children
        li = tuple(sorted(positions[x] for x in descendants[left]))
        ri = tuple(sorted(positions[x] for x in descendants[right]))
        eigen = eigendecompose_covariance(whitened[list(li + ri)], compute_eigenvectors=True)
        if eigen is None:
            raise ValueError("Discovery node has no positive PCA directions.")
        projection, eigenvalues = build_pca_projection(
            eigen,
            projection_dimension=values.shape[1],
            n_features_total=values.shape[1],
        )
        if not np.allclose(projection @ projection.T, np.eye(len(projection)), atol=1e-8):
            raise ValueError("Discovery PCA directions are not orthonormal.")
        nodes.append(
            IndependentGaussianNode(
                parent,
                left,
                right,
                li,
                ri,
                _readonly(projection),
                _readonly(eigenvalues),
            )
        )
    return IndependentGaussianTree(
        nx.freeze(tree),
        tuple(discovery.index),
        tuple(discovery.columns),
        ids,
        _digest(values),
        _readonly(cholesky),
        tuple(nodes),
        energy_fraction,
        frozenset(tree.edges),
        tuple((u, v, float(attrs["branch_length"])) for u, v, attrs in tree.edges(data=True)),
        diffusion_coordinates,
        json.dumps(geometry_metadata, sort_keys=True),
        json.dumps(branch_metadata, sort_keys=True),
    )


def test_independent_gaussian_tree(
    model: IndependentGaussianTree,
    inference: pd.DataFrame,
    *,
    observation_ids,
    alpha: float = 0.01,
    edge_alpha: float = 0.001,
    passthrough: bool = True,
) -> IndependentGaussianResult:
    """Replace nominal/empirical tails on a frozen independent Gaussian tree.

    Adaptive ranks use the conditional energy law. Fixed ranks retain their
    exact chi-square reference. Both child-parent contrasts at a binary parent
    standardize to the same sibling contrast (up to sign) under known common
    covariance; therefore their raw p-values coincide. Edge and sibling multiple
    testing still run through their separate existing correction procedures.
    """
    values = _values(inference)
    if (
        tuple(inference.index) != model.leaf_labels
        or tuple(inference.columns) != model.feature_labels
    ):
        raise ValueError("Inference leaf order and feature labels must match discovery exactly.")
    ids = _observations(observation_ids, len(values))
    if ids & model.discovery_observation_ids:
        raise ValueError("Discovery and inference observation IDs must be disjoint.")
    digest = _digest(values)
    if digest == model.discovery_sha256:
        raise ValueError("Inference reuses the discovery measurements.")
    if frozenset(model.tree.edges) != model.tree_edges:
        raise ValueError("Discovery topology changed after fitting.")
    if any(
        model.tree.edges[u, v].get("branch_length") != length
        for u, v, length in model.tree_branch_lengths
    ):
        raise ValueError("Discovery branch lengths changed after fitting.")
    if model.diffusion_coordinates is not None:
        if (
            _digest(model.diffusion_coordinates)
            != json.loads(model.geometry_metadata_json)["coordinates_sha256"]
        ):
            raise ValueError("Discovery diffusion coordinates changed after fitting.")
    if any(not np.isfinite(a) or not 0 < a <= 1 for a in (alpha, edge_alpha)):
        raise ValueError("Gate alpha values must be finite and in (0, 1].")
    whitened = solve_triangular(model.covariance_cholesky, values.T, lower=True).T
    rows = []
    for node in model.nodes:
        left, right = whitened[list(node.left_rows)], whitened[list(node.right_rows)]
        z = (left.mean(0) - right.mean(0)) / np.sqrt(1 / len(left) + 1 / len(right))
        nominal = run_projected_wald_kernel(
            z,
            spectral_k=len(node.projection),
            pca_projection=node.projection,
            pca_eigenvalues=node.eigenvalues,
            adaptive_dimension_energy_fraction=model.energy_fraction,
        )
        k = nominal.projection_dimension
        previous = k - 1
        ties = tied_eigenvalue_neighbors(node.eigenvalues)
        while previous > 0 and ties[previous - 1]:
            previous -= 1
        p_value = (
            nominal.p_value
            if model.energy_fraction is None
            else energy_selected_gaussian_sf(
                nominal.statistic,
                candidate_dimension=len(node.projection),
                projection_dimension=k,
                previous_group_dimension=previous,
                energy_fraction=model.energy_fraction,
            )
        )
        rows.append(
            dict(
                parent=node.parent,
                left=node.left,
                right=node.right,
                left_n=len(left),
                right_n=len(right),
                statistic=nominal.statistic,
                projection_dimension=k,
                candidate_dimension=len(node.projection),
                previous_group_dimension=previous,
                energy_fraction=model.energy_fraction,
                nominal_p_value=nominal.p_value,
                p_value=p_value,
                reference_law="fixed_gaussian_chi_square"
                if model.energy_fraction is None
                else "energy_selected_gaussian_conditional",
            )
        )
    tests = pd.DataFrame(rows).set_index("parent", drop=False)
    annotations, decomposition = _annotate_tests(
        model, tests, alpha=alpha, edge_alpha=edge_alpha, passthrough=passthrough
    )
    return IndependentGaussianResult(tests, annotations, decomposition, digest)


def _annotate_tests(model, tests, *, alpha, edge_alpha, passthrough):
    """Shared correction/traversal stage; benchmark rescoring holds tests fixed."""
    edge_children = [child for row in tests.itertuples() for child in (row.left, row.right)]
    edge_p = [row.p_value for row in tests.itertuples() for _ in range(2)]
    edge = apply_tree_bh_correction(model.tree, np.array(edge_p), edge_children, alpha=edge_alpha)
    annotations = init_sibling_annotation_df(
        pd.DataFrame(
            {"Child_Parent_Divergence_Significant": False},
            index=list(model.tree.nodes),
        )
    )
    annotations.loc[edge_children, "Child_Parent_Divergence_Significant"] = (
        edge.child_parent_edge_null_rejected_by_tree_bh
    )
    annotations.loc[edge_children, "Child_Parent_Divergence_P_Value"] = edge_p
    annotations.loc[edge_children, "Child_Parent_Divergence_P_Value_Corrected"] = (
        edge.child_parent_edge_corrected_p_values_by_tree_bh
    )
    annotations = apply_traversal_aligned_sibling_bh_results(
        model.tree,
        annotations,
        list(tests.index),
        [(r.statistic, r.projection_dimension, r.p_value) for r in tests.itertuples()],
        alpha=alpha,
    )
    # Existing transport column carries retained rank; it is not a chi-square
    # reference df when the reference is the energy-selected mixture.
    annotations.loc[tests.index, "Sibling_Projection_Dimension"] = tests.projection_dimension
    annotations.loc[tests.index, "Sibling_Test_Method"] = tests.reference_law
    annotations.loc[tests.index, "Sibling_Gate_P_Value_Calibration"] = tests.reference_law
    annotations.loc[tests.index, "Sibling_Gate_P_Value_Role"] = "active_traversal_sibling_gate"
    annotations.loc[tests.index, "Sibling_Degrees_of_Freedom_Semantics"] = (
        "retained_rank_not_fitted_reference_df"
    )
    decomposition = TreeDecomposition(
        tree=model.tree,
        annotations_df=annotations,
        passthrough=passthrough,
        trace_level="full",
    ).decompose_tree()
    return annotations, decomposition


# Avoid accidental pytest collection when imported into a test module.
test_independent_gaussian_tree.__test__ = False

__all__ = [
    "IndependentGaussianTree",
    "IndependentGaussianResult",
    "IndependentGaussianDiffusion",
    "fit_independent_gaussian_tree",
    "test_independent_gaussian_tree",
]
