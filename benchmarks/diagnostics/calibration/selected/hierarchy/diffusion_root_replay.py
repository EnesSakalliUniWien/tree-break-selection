"""Serial observable Gaussian radial replay of actual diffusion; no cell certificates.

Every evaluation refits neighbors, KDE, epsilon, spectrum and average linkage.
Instrumentation delegates unchanged to native implementations. It temporarily
observes module functions, so run this diagnostic serially, outside concurrent
production geometry calls. Point margins are not uniform interval bounds.
"""

import hashlib
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from pydiffmap import kernel
from scipy.spatial.distance import squareform

from tree_break_selection.eigenvalue_groups import tied_eigenvalue_neighbors
from tree_break_selection.space_separation import diffusion
from tree_break_selection.tree.construction.build import _canonical_sample_order, build_tree
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns


@dataclass(frozen=True)
class DiffusionReplayPoint:
    radius: float
    status: str
    signature: tuple
    topology_signature: tuple
    margins: np.ndarray
    comparisons: tuple
    state: dict
    error: str | None = None


def _topology(signature, n):
    groups = {i: (i,) for i in range(n)}
    clades = []
    for step, (a, b) in enumerate(signature):
        groups[n + step] = tuple(sorted(groups[a] + groups[b]))
        clades.append(groups[n + step])
    return tuple(sorted(clades))


def _geometry(values, labels):
    state = {}
    original_kernel, original_eigh = kernel.Kernel, diffusion.eigh
    original_epsilon = diffusion.resolve_adaptive_epsilon

    def factory(*args, **kwargs):
        fitted = original_kernel(*args, **kwargs)
        state["kernel"] = fitted
        return fitted

    def epsilon(fitted, *args, **kwargs):
        # Snapshot before Kernel.compute mutates its sparse distances into weights.
        state["scaled_distances"] = fitted.scaled_dists.copy()
        return original_epsilon(fitted, *args, **kwargs)

    def eigh(*args, **kwargs):
        result = original_eigh(*args, **kwargs)
        state["eigenvalues"] = result[0][::-1].copy()
        return result

    data = pd.DataFrame(values, index=labels)
    space = continuous_feature_space_from_columns(data.columns)
    with (
        patch("pydiffmap.kernel.Kernel", factory),
        patch.object(diffusion, "resolve_adaptive_epsilon", epsilon),
        patch.object(diffusion, "eigh", eigh),
    ):
        geometry = diffusion.feature_family_diffusion_geometry(data, space)
    fitted = state.pop("kernel")
    neighbor_distances, neighbor_indices = fitted.neigh.kneighbors(values)
    scaled_snapshot = state["scaled_distances"]
    diagonal = scaled_snapshot.diagonal()
    all_positive = scaled_snapshot.data[scaled_snapshot.data > 0.0]
    off_diagonal = scaled_snapshot.copy()
    off_diagonal.setdiag(0.0)
    off_diagonal.eliminate_zeros()
    state.update(
        positive_self_scaled_distance_count=int(np.count_nonzero(diagonal > 0.0)),
        zero_self_scaled_distance_count=int(np.count_nonzero(diagonal == 0.0)),
        positive_scaled_distance_count=int(len(all_positive)),
        actual_positive_scaled_sq_median=float(np.median(all_positive**2)),
        off_diagonal_positive_scaled_sq_median=float(
            np.median(off_diagonal.data[off_diagonal.data > 0.0] ** 2)
        ),
    )
    state.update(
        neighbors=neighbor_indices.copy(),
        neighbor_distances=neighbor_distances.copy(),
        bandwidths=np.array(fitted.bandwidths, copy=True),
        kde_dimension=float(fitted.kde.d),
        kde_epsilon=float(fitted.kde.epsilon_fitted),
        kde_bgh_index=int(np.argmin(abs(2.0 ** np.arange(-40, 41) - fitted.kde.epsilon_fitted))),
        outer_epsilon=float(geometry.metadata["epsilon"]),
        retained_dimension=geometry.coordinates.shape[1],
        metadata=dict(geometry.metadata),
    )
    if "eigenvalues" in state:
        eigenvalues = state["eigenvalues"]
        start, end = diffusion._leading_diffusion_spectrum(
            eigenvalues, min(50, len(values) - 2, values.shape[1])
        )
        state["retained_indices"] = (start, end)
        state["eigenvalue_ties"] = tied_eigenvalue_neighbors(eigenvalues)
    else:
        state["spectral_observation"] = "full spectrum unavailable on native sparse solver route"
    built = build_tree(data, geometry.distance_condensed, linkage_method="average")
    signature = tuple(tuple(sorted(map(int, row[:2]))) for row in built.linkage_matrix)
    state["distance_condensed"] = geometry.distance_condensed.copy()
    return signature, state


def _margins(distances, observed_signature):
    n = distances.shape[0]
    active = set(range(n))
    counts = {i: 1 for i in active}
    scores = {(a, b): float(distances[a, b]) for a, b in combinations(range(n), 2)}
    margins = []
    comparisons = []
    for step, chosen in enumerate(observed_signature):
        selected_score = scores[chosen]
        for competitor in combinations(sorted(active), 2):
            if competitor == chosen:
                continue
            margins.append(scores[competitor] - selected_score)
            comparisons.append((step, chosen, competitor))
        a, b = chosen
        merged = n + step
        counts[merged] = counts[a] + counts[b]
        for other in active - {a, b}:
            scores[tuple(sorted((merged, other)))] = (
                counts[a] * scores[tuple(sorted((a, other)))]
                + counts[b] * scores[tuple(sorted((b, other)))]
            ) / counts[merged]
        active.difference_update(chosen)
        active.add(merged)
    return np.array(margins), tuple(comparisons)


@dataclass(frozen=True)
class DiffusionRootPath:
    values: np.ndarray
    labels: tuple
    offset: np.ndarray
    slope: np.ndarray
    e: np.ndarray
    u: np.ndarray
    sigma: float
    dimension: int
    observed_radius: float
    observed_signature: tuple
    source_hashes: dict
    tie_tolerance: float = 1e-12

    def evaluate(self, radius):
        """Observe a point; ties and geometry failures remain explicit unresolved states."""
        radius = float(radius)
        if not np.isfinite(radius) or radius < 0:
            raise ValueError("radius must be finite and nonnegative")
        values = self.offset + radius * self.slope
        try:
            signature, state = _geometry(values, self.labels)
        except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
            return DiffusionReplayPoint(
                radius,
                "unsupported",
                (),
                (),
                np.array([]),
                (),
                {},
                f"{type(error).__name__}: {error}",
            )
        distances = squareform(state["distance_condensed"])
        margins, comparisons = _margins(distances, self.observed_signature)
        tolerance = self.tie_tolerance * max(1.0, float(distances.max()))
        state["merge_tie_tolerance"] = tolerance
        state["near_tie_comparisons"] = np.flatnonzero(abs(margins) <= tolerance)
        status = "selected" if signature == self.observed_signature else "nonselected"
        if len(state["near_tie_comparisons"]):
            status = "tie"
        elif (margins > 0).all() != (status == "selected"):
            status = "unsupported"
            state["inconsistency"] = "native signature and recomputed merge comparisons disagree"
        return DiffusionReplayPoint(
            radius,
            status,
            signature,
            _topology(signature, len(values)),
            margins,
            comparisons,
            state,
        )


def make_diffusion_root_path(values, *, labels=None, sigma=1.0, tie_tolerance=1e-12):
    """Canonicalize labels once, then build the full-feature known-I radial path.

    This diagnostic uses the native constructor's canonical-label ordering before
    geometry fitting, and preserves that order at every radius. Input is copied.
    Observed unsupported/tied data cannot define a supported replay path here.
    """
    if isinstance(values, pd.DataFrame):
        if labels is not None:
            raise ValueError("supply labels through the DataFrame index only")
        labels = values.index.tolist()
        values = values.to_numpy(dtype=float)
    x = np.array(values, dtype=float, copy=True)
    if x.ndim != 2 or len(x) < 2 or x.shape[1] < 1 or not np.isfinite(x).all():
        raise ValueError("values must be finite n-by-p with n>=2,p>=1")
    if not np.isscalar(sigma) or not np.isfinite(sigma) or sigma <= 0:
        raise ValueError("sigma must be known finite positive scalar")
    if not np.isfinite(tie_tolerance) or tie_tolerance < 0:
        raise ValueError("tie_tolerance must be finite and nonnegative")
    labels = list(labels) if labels is not None else [f"S{i}" for i in range(len(x))]
    if len(labels) != len(x) or not pd.Index(labels).is_unique:
        raise ValueError("labels must be unique and match row count")
    order = _canonical_sample_order(labels)
    x = x[order]
    labels = tuple(labels[i] for i in order)
    try:
        signature, _ = _geometry(x, labels)
    except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
        raise ValueError(f"observed diffusion geometry unsupported: {error}") from error
    groups = {i: [i] for i in range(len(x))}
    for step, (a, b) in enumerate(signature):
        groups[len(x) + step] = groups[a] + groups[b]
    a, b = signature[-1]
    eta = np.zeros(len(x))
    eta[groups[a]] = 1 / len(groups[a])
    eta[groups[b]] = -1 / len(groups[b])
    e = eta / np.linalg.norm(eta)
    contrast = e @ x
    radius = np.linalg.norm(contrast) / sigma
    if not radius > 0:
        raise ValueError("observed contrast direction undefined")
    u = contrast / np.linalg.norm(contrast)
    offset = x - np.outer(e, contrast)
    slope = sigma * np.outer(e, u)
    sources = [
        Path(diffusion.__file__),
        Path(kernel.__file__),
        Path(__file__),
        Path(build_tree.__code__.co_filename),
        Path(tied_eigenvalue_neighbors.__code__.co_filename),
    ]
    hashes = {str(p.resolve()): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    path = DiffusionRootPath(
        x,
        labels,
        offset,
        slope,
        e,
        u,
        float(sigma),
        x.shape[1],
        float(radius),
        signature,
        hashes,
        float(tie_tolerance),
    )
    observed = path.evaluate(radius)
    if observed.status != "selected":
        raise ValueError(f"observed replay status is {observed.status}: {observed.error}")
    return path
