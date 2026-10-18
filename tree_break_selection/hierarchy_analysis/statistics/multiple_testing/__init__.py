"""Multiple testing correction utilities for statistical hypothesis testing.

This package provides multiple-testing utilities used by the Tree-BH
edge-divergence pipeline.

Modules
-------
base
    Core Benjamini-Hochberg FDR correction
tree_bh/
    Hierarchical Tree-BH package (Bogomolov et al. 2021)
"""

from .base import benjamini_hochberg_correction
from .tree_bh import apply_tree_bh_correction

__all__ = [
    # Core function
    "benjamini_hochberg_correction",
    # Tree-BH correction
    "apply_tree_bh_correction",
]
