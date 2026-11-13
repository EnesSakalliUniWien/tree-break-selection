"""The table builder's search finds the best flat cut, which is the best of all splits."""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

import build_best_flat_cut_reference as builder  # noqa: E402

from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.best_flat_cut_reference import (  # noqa: E402
    split_rho,
)


def _best_of_all_splits(x: np.ndarray) -> float:
    n = len(x)
    best = 0.0
    for bits in itertools.product([False, True], repeat=n - 1):
        in_left = np.array((True, *bits))
        if in_left.all():
            continue
        best = max(best, split_rho(x, in_left))
    return best


def test_search_equals_the_maximum_over_all_splits_on_small_datasets() -> None:
    rng = np.random.default_rng(4)
    for k in (2, 3):
        for _ in range(10):
            x = rng.standard_normal((10, k))
            assert np.isclose(builder.best_flat_cut_rho(x), _best_of_all_splits(x), rtol=1e-9)
