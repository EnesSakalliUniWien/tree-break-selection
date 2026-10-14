from __future__ import annotations

import sys

import numpy as np
import pytest
from benchmarks.shared.plots.manifold import create_manifold_alignment_plot


def test_manifold_plot_fails_without_umap_instead_of_relabelling_pca(monkeypatch) -> None:
    monkeypatch.setitem(sys.modules, "umap", None)
    data = np.random.default_rng(0).normal(size=(12, 4))
    with pytest.raises(ImportError):
        create_manifold_alignment_plot(data, np.zeros(12, dtype=int), 1, {})
