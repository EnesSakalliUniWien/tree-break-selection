import importlib
import importlib.util

import networkx as nx
import numpy as np
import pandas as pd
import pytest


def _expand_leaves_into_rows(tree: nx.DiGraph, _annotations: pd.DataFrame, leaf_data: pd.DataFrame):
    """Hang four identical data rows below every leaf and rebuild counts and means.

    The edge test fits its basis on one half of a parent's rows and tests the other
    half, so a parent whose children are single rows has no edge test.
    """
    tree = tree.copy()
    rows = {}
    for leaf in leaf_data.index:
        row = leaf_data.loc[leaf].to_numpy(dtype=np.float64)
        for copy in range(4):
            rows[f"{leaf}{copy}"] = row
            tree.add_edge(leaf, f"{leaf}{copy}", branch_length=0.05)
            tree.nodes[f"{leaf}{copy}"].update(is_leaf=True, label=f"{leaf}{copy}")
        tree.nodes[leaf]["is_leaf"] = False
    for node in reversed(list(nx.topological_sort(tree))):
        children = list(tree.successors(node))
        if not children:
            tree.nodes[node].update(leaf_count=1, distribution=rows[node])
            continue
        counts = [tree.nodes[child]["leaf_count"] for child in children]
        tree.nodes[node]["leaf_count"] = sum(counts)
        tree.nodes[node]["distribution"] = sum(
            count * tree.nodes[child]["distribution"] for count, child in zip(counts, children)
        ) / sum(counts)
    annotations = pd.DataFrame(
        {"leaf_count": {node: tree.nodes[node]["leaf_count"] for node in tree.nodes}}
    )
    return (
        tree,
        annotations,
        pd.DataFrame.from_dict(rows, orient="index", columns=leaf_data.columns),
    )


@pytest.fixture
def expand_leaves_into_rows():
    """Return a helper that turns a (tree, annotations, leaf_data) fixture into row-level data."""
    return _expand_leaves_into_rows


@pytest.fixture(scope="module")
def require_optional_dependencies():
    """Import optional modules, skipping only when they are not installed.

    An installed module whose import fails is a broken optional environment and
    must remain visible as a test failure instead of being silently skipped.
    """

    def require(*module_names: str) -> tuple[object, ...]:
        missing = [name for name in module_names if importlib.util.find_spec(name) is None]
        if missing:
            pytest.skip(f"Missing optional dependencies: {', '.join(missing)}")
        return tuple(importlib.import_module(name) for name in module_names)

    return require
