"""Opt-in sibling p-value: a split's ``rho`` against a Gaussian reference for the best flat cut.

For a split of a subtree's n cells into left and right, inside the subtree's k invariant directions
(``invariant_directions``), ``rho = pi(1-pi) D^2`` with the pooled within covariance (divisor n - 2),
pi the left share and D the Mahalanobis distance of the two sides' centroids.

Exact:
- For one Gaussian population, ``rho`` of the best flat cut in k fixed directions does not change under
  shifts and invertible linear maps, so its distribution depends only on (n, k).
- The best flat cut is the best of all 2^(n-1) - 1 splits (checked by enumeration for n = 10 to 14), so any
  split's own ``rho`` is at most the best flat cut's ``rho``. Comparing a split with the null of the best
  flat cut therefore needs no count of splits.

Only measured (wiki page ``gaussian-best-flat-cut-reference``):
- The stored null sample comes from a search (re-sorts from the top 50 directions plus single-point moves),
  not from exact enumeration; it matches enumeration for n up to 40 and a 300,000-direction search up to
  k = 6.
- n between grid values is rounded down to the largest grid n, which errs on the strict side because the null
  falls with n; this is checked on the table, not proved.
- k comes from ``invariant_directions`` on the same data. Measured: no excess rejections at n = 400; a small
  excess at n = 100 with 50 or 80 features (152/12,000 at the 99% point).
- For k above 10 the search is known not to find the best flat cut reliably (checked at k = 15 and 20: its
  99% points are lower bounds). The stored null for k above 10 is therefore the null of the searched value,
  not of the best flat cut, and for an arbitrary split (the tree's split) the p-value is not guaranteed
  conservative there. Measured under one Gaussian population at k = 11, 13, 15 (wiki page
  ``best-flat-cut-reference-k15``): the root split of a Ward tree exceeded the 99% point in 0 of 12,000
  datasets; the best of stronger searches exceeded it in 25 of 2,700 (0.93%), 6 of 300 at n = 800, k = 15.
- Limits: k from 1 to ``MAX_DIRECTIONS`` = 15 and n from ``MIN_CELLS`` = 15, and for k above 10 n from
  ``MIN_CELLS_ABOVE_10_DIRECTIONS`` = 50 (those cells are not built); outside them the p-value is 1.
- A group below ``MIN_CELLS`` can never be tested, so with the opt-in the gate pipeline lets a small peel
  (smaller child below ``MIN_CELLS`` leaves, larger child at least ``MIN_CELLS``;
  ``tree.topology.peel_continuations``) pass the edge gate like a singleton peel. Measured (wiki page
  ``k-split-direction-draft``): one Gaussian population unchanged in 880 of 880 fits; CHOIR cells, default
  tree, ARI 0.0275 to 0.6210; Markov/Ward tree 0.8219 to 0.8214 with 116 to 110 cell lines recovered.

The p-value is a Monte Carlo p-value, (1 + #{null >= rho}) / (B + 1), so it cannot go below 1 / (B + 1).
Only the largest ``ceil(0.1 * B)`` null values per (n, k) are stored; a ``rho`` below the smallest stored
value gets p = 1, which never matters for a Benjamini-Hochberg rejection at alpha below 0.1.

Build the table with ``scripts/build_best_flat_cut_reference.py``.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

from .split_count_rule_out import _leaf_rows, invariant_directions

BEST_FLAT_CUT_REFERENCE_LABEL = "best_flat_cut_reference"
MIN_CELLS = 15
MIN_CELLS_ABOVE_10_DIRECTIONS = 50
MAX_DIRECTIONS = 15
TABLE_PATH = Path(__file__).with_name("best_flat_cut_reference_table.npz")


@dataclass(frozen=True)
class BestFlatCutReference:
    """Reference p-value of one split and why it was tested or failed closed."""

    p_value: float
    reason: str  # tested, no_invariant_directions, too_many_invariant_directions or too_few_cells
    invariant_directions: int
    rho: float  # nan unless tested
    reference_n: int  # largest grid n <= n; 0 unless tested


def whiten(points: np.ndarray) -> np.ndarray:
    """Centre ``points`` and whiten them by the Cholesky factor of their covariance (divisor n)."""
    n = len(points)
    centred = points - points.mean(0)
    return centred @ np.linalg.inv(np.linalg.cholesky(centred.T @ centred / n)).T


def to_rho(r2: float, n: int) -> float:
    """Return ``rho`` from the between/total ratio ``r2`` of a split of n cells; inf when r2 >= 1."""
    return float("inf") if r2 >= 1.0 else float((n - 2) / n * r2 / (1.0 - r2))


def split_rho(coordinates: np.ndarray, in_left: np.ndarray) -> float:
    """Return ``rho`` of a split of the cells' ``coordinates`` (n x k) into ``in_left`` and the rest."""
    in_left = np.asarray(in_left, dtype=bool)
    if in_left.all() or not in_left.any():
        raise ValueError("Best-flat-cut reference needs cells on both sides of the split.")
    n, size = len(coordinates), int(in_left.sum())
    whitened = whiten(np.asarray(coordinates, dtype=np.float64))
    r2 = float((whitened[in_left].sum(0) ** 2).sum() / (size * (n - size)))
    return to_rho(r2, n)


@lru_cache(maxsize=1)
def load_table() -> dict[str, np.ndarray]:
    """Load the stored upper null sample: ``upper[i, j]`` is sorted for (n_grid[i], k_grid[j])."""
    with np.load(TABLE_PATH) as stored:
        return {name: stored[name] for name in stored.files}


def reference_p_value(rho: float, n: int, k: int) -> tuple[float, int]:
    """Return the Monte Carlo p-value of ``rho`` at (n, k) and the grid n used; n >= 15, 1 <= k <= 15.

    For k > 10, n must be at least ``MIN_CELLS_ABOVE_10_DIRECTIONS``; the cells below are NaN.
    """
    table = load_table()
    n_grid = table["n_grid"]
    i = int(np.searchsorted(n_grid, n, side="right")) - 1
    j = int(np.flatnonzero(table["k_grid"] == k)[0])
    upper = table["upper"][i, j].astype(np.float64)
    if np.isnan(upper).any():
        raise ValueError(f"No reference for n = {n}, k = {k}.")
    if rho < upper[0]:
        return 1.0, int(n_grid[i])
    count = len(upper) - int(np.searchsorted(upper, rho, side="left"))
    return (1 + count) / (int(table["B"]) + 1), int(n_grid[i])


def best_flat_cut_reference(points: np.ndarray, in_left: np.ndarray) -> BestFlatCutReference:
    """Return the reference p-value of a split of ``points`` into ``in_left`` and the rest."""
    points = np.asarray(points, dtype=np.float64)
    in_left = np.asarray(in_left, dtype=bool)
    if in_left.all() or not in_left.any():
        raise ValueError("Best-flat-cut reference needs cells on both sides of the split.")
    n = len(points)
    directions = invariant_directions(points)
    k = len(directions)
    if k == 0:
        reason = "no_invariant_directions"
    elif k > MAX_DIRECTIONS:
        reason = "too_many_invariant_directions"
    elif n < MIN_CELLS or (k > 10 and n < MIN_CELLS_ABOVE_10_DIRECTIONS):
        reason = "too_few_cells"
    else:
        rho = split_rho(points @ directions.T, in_left)
        p_value, reference_n = reference_p_value(rho, n, k)
        return BestFlatCutReference(p_value, "tested", k, rho, reference_n)
    return BestFlatCutReference(1.0, reason, k, float("nan"), 0)


def apply_best_flat_cut_reference(
    tree: nx.DiGraph,
    leaf_data: pd.DataFrame,
    records: list,
    tested_parent_ids: list,
    summaries: list[tuple[float, float, float]],
    method_labels: list[str],
    decisions: list,
) -> tuple[list[tuple[float, float, float]], list[str], list]:
    """Swap each tested split's p-value for its reference p-value; zero-dimensional records stay.

    The summary is (rho, k, p-value). Downstream requires a finite statistic, so a split whose ``rho`` is
    not finite (failed closed, or perfectly separated) keeps its projected-Wald statistic in that slot.
    """
    row_of_label = {label: row for row, label in enumerate(leaf_data.index)}
    matrix = leaf_data.to_numpy(dtype=np.float64)
    by_parent = {record.parent: record for record in records}
    new_summaries, new_labels, new_decisions = [], [], []
    for parent, summary, label, decision in zip(
        tested_parent_ids, summaries, method_labels, decisions, strict=True
    ):
        if decision is None:
            new_summaries.append(summary)
            new_labels.append(label)
            new_decisions.append(decision)
            continue
        record = by_parent[parent]
        left = _leaf_rows(tree, record.left, row_of_label)
        right = _leaf_rows(tree, record.right, row_of_label)
        result = best_flat_cut_reference(
            matrix[left + right],
            np.r_[np.ones(len(left), dtype=bool), np.zeros(len(right), dtype=bool)],
        )
        statistic = (
            result.rho
            if np.isfinite(result.rho)
            else float(record.stat) / float(record.reference_scale)
        )
        new_summaries.append((statistic, float(result.invariant_directions), result.p_value))
        new_labels.append(f"{BEST_FLAT_CUT_REFERENCE_LABEL}:{result.reason}")
        new_decisions.append(replace(decision, c_hat=1.0, p_value=result.p_value))
    return new_summaries, new_labels, new_decisions


__all__ = [
    "BEST_FLAT_CUT_REFERENCE_LABEL",
    "MAX_DIRECTIONS",
    "MIN_CELLS",
    "MIN_CELLS_ABOVE_10_DIRECTIONS",
    "BestFlatCutReference",
    "apply_best_flat_cut_reference",
    "best_flat_cut_reference",
    "load_table",
    "reference_p_value",
    "split_rho",
    "to_rho",
    "whiten",
]
