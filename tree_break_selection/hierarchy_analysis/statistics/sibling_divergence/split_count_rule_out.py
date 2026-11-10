"""Opt-in sibling p-value: count every possible split, then rule out splits that are not real.

A tree chooses each split among the 2^(n-1) - 1 two-group splits of a subtree's n cells, so the
p-value is the union bound min(1, (2^(n-1) - 1) * chi2.sf(T, df)) over all of them. A split is also
ruled out (p = 1) unless

1. a flat cut in the subtree's invariant directions separates its two sides. Invariant directions are
   eigenvectors of the subtree's cell covariance whose eigenvalue lies above its noise band: noise
   level = median eigenvalue / median of the Marchenko-Pastur law with ratio p/n, band edge =
   level * (1 + sqrt(p/n))^2. A subtree of n <= p + 1 cells has none; and
2. at most ``MAX_EQUIDISTANT_SHARE`` of its cells are about equally far from both sides' centroids in
   the test metric (cells whitened by the pooled within-sides covariance): |margin| <
   ``EQUIDISTANT_MARGIN``, margin = (distance to the other centroid - distance to its own) / distance
   between the centroids.

Measured on simulated and real data with 30 test coordinates (wiki pages split-count-exclusion and
split-count-exclusion-real). The exclusions look at the observed split, so error control is measured, not
guaranteed. Part of the protection of small subtrees comes from rule 1 giving no invariant directions at
n <= p + 1; with fewer coordinates (for example 10) small noise subtrees keep invariant directions and their
splits can pass.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import networkx as nx
import numpy as np
import pandas as pd
from scipy.optimize import linprog
from scipy.stats import chi2

MAX_EQUIDISTANT_SHARE = 0.05
EQUIDISTANT_MARGIN = 0.1
SPLIT_COUNT_RULE_OUT_LABEL = "split_count_rule_out"


@dataclass(frozen=True)
class SplitCountRuleOut:
    """Rule-out p-value of one split and why it was kept or ruled out."""

    p_value: float
    reason: str  # kept, no_invariant_directions, not_separable or equidistant_cells
    invariant_directions: int
    equidistant_share: float


def log_split_count(n: int) -> float:
    """Return log(2^(n-1) - 1), the number of two-group splits of n cells (0 for n <= 2)."""
    if n <= 2:
        return 0.0
    bits = (n - 1) * np.log(2.0)
    return float(bits + np.log(-np.expm1(-bits)))


def _marchenko_pastur_median(ratio: float) -> float:
    lower, upper = (1 - np.sqrt(ratio)) ** 2, (1 + np.sqrt(ratio)) ** 2
    grid = np.linspace(lower, upper, 20001)
    density = np.sqrt(np.maximum((upper - grid) * (grid - lower), 0.0)) / (2 * np.pi * ratio * grid)
    cumulative = np.cumsum(density)
    return float(grid[np.searchsorted(cumulative / cumulative[-1], 0.5)])


def invariant_directions(points: np.ndarray) -> np.ndarray:
    """Return the eigenvectors (rows) of the cells' covariance above the noise band of its bulk."""
    n, p = points.shape
    if n <= p + 1:
        return np.zeros((0, p))
    values, vectors = np.linalg.eigh(np.cov(points.T))
    values, vectors = values[::-1], vectors[:, ::-1]
    ratio = p / n
    edge = np.median(values) / _marchenko_pastur_median(ratio) * (1 + np.sqrt(ratio)) ** 2
    return vectors[:, values > edge].T


def is_separable(points: np.ndarray, in_left: np.ndarray) -> bool:
    """Is there a flat cut w.x + b with every left cell on one side and every right cell on the other?"""
    sign = np.where(in_left, 1.0, -1.0)
    result = linprog(
        np.zeros(points.shape[1] + 1),
        A_ub=-np.column_stack([sign[:, None] * points, sign]),
        b_ub=-np.ones(len(sign)),
        bounds=(None, None),
        method="highs",
    )
    if result.status not in (0, 2):
        raise RuntimeError(f"Separability linear program failed: {result.message}")
    return result.status == 0


def equidistant_share(points: np.ndarray, in_left: np.ndarray) -> float:
    """Share of cells about equally far from both sides' centroids in the whitened test metric."""
    left, right = points[in_left], points[~in_left]
    centered = np.vstack([left - left.mean(0), right - right.mean(0)])
    values, vectors = np.linalg.eigh(centered.T @ centered / max(len(points) - 2, 1))
    whitened = points @ vectors / np.sqrt(np.maximum(values, values.max() * 1e-12))
    mu_left, mu_right = whitened[in_left].mean(0), whitened[~in_left].mean(0)
    to_left = np.linalg.norm(whitened - mu_left, axis=1)
    to_right = np.linalg.norm(whitened - mu_right, axis=1)
    margin = np.where(in_left, to_right - to_left, to_left - to_right) / np.linalg.norm(
        mu_left - mu_right
    )
    return float(np.mean(np.abs(margin) < EQUIDISTANT_MARGIN))


def split_count_rule_out(
    statistic: float,
    degrees_of_freedom: float,
    points: np.ndarray,
    in_left: np.ndarray,
) -> SplitCountRuleOut:
    """Return the rule-out p-value of a split of ``points`` into ``in_left`` and the rest."""
    points = np.asarray(points, dtype=np.float64)
    in_left = np.asarray(in_left, dtype=bool)
    if in_left.all() or not in_left.any():
        raise ValueError("Split-count rule-out needs cells on both sides of the split.")
    p_value = float(
        np.exp(min(0.0, log_split_count(len(points)) + chi2.logsf(statistic, degrees_of_freedom)))
    )
    directions = invariant_directions(points)
    share = float(
        "nan"
    )  # measured only for separable splits; it needs a full-rank p x p covariance
    if len(directions) == 0:
        reason = "no_invariant_directions"
    elif not is_separable(points @ directions.T, in_left):
        reason = "not_separable"
    else:
        share = equidistant_share(points, in_left)
        reason = "equidistant_cells" if share > MAX_EQUIDISTANT_SHARE else "kept"
    return SplitCountRuleOut(
        p_value=p_value if reason == "kept" else 1.0,
        reason=reason,
        invariant_directions=len(directions),
        equidistant_share=share,
    )


def _leaf_rows(tree: nx.DiGraph, node: object, row_of_label: dict[object, int]) -> list[int]:
    below = nx.descendants(tree, node) | {node}
    return [row_of_label[tree.nodes[u]["label"]] for u in below if tree.out_degree(u) == 0]


def apply_split_count_rule_out(
    tree: nx.DiGraph,
    leaf_data: pd.DataFrame,
    records: list,
    tested_parent_ids: list,
    summaries: list[tuple[float, float, float]],
    method_labels: list[str],
    decisions: list,
) -> tuple[list[tuple[float, float, float]], list[str], list]:
    """Swap each tested split's p-value for its rule-out p-value; zero-dimensional records stay."""
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
        statistic = float(record.stat) / float(record.reference_scale)
        result = split_count_rule_out(
            statistic,
            float(record.degrees_of_freedom),
            matrix[left + right],
            np.r_[np.ones(len(left), dtype=bool), np.zeros(len(right), dtype=bool)],
        )
        new_summaries.append((statistic, float(record.degrees_of_freedom), result.p_value))
        new_labels.append(f"{SPLIT_COUNT_RULE_OUT_LABEL}:{result.reason}")
        new_decisions.append(replace(decision, c_hat=1.0, p_value=result.p_value))
    return new_summaries, new_labels, new_decisions


__all__ = [
    "EQUIDISTANT_MARGIN",
    "MAX_EQUIDISTANT_SHARE",
    "SPLIT_COUNT_RULE_OUT_LABEL",
    "SplitCountRuleOut",
    "apply_split_count_rule_out",
    "equidistant_share",
    "invariant_directions",
    "is_separable",
    "log_split_count",
    "split_count_rule_out",
]
