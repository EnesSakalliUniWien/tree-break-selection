"""Contracts for converting explicit cluster boundaries to assignments."""

from __future__ import annotations

import pytest
from tree_break_selection.hierarchy_analysis.cluster_assignments import (
    ClusterBoundary,
    build_cluster_assignments,
    build_sample_cluster_assignments,
)


def test_cluster_assignments_use_explicit_boundary_root() -> None:
    assignments = build_cluster_assignments(
        [
            ClusterBoundary(root_node="left_subtree", leaves=frozenset({"L1", "L2"})),
            ClusterBoundary(root_node="right_leaf", leaves=frozenset({"R"})),
        ]
    )

    assert assignments == {
        0: {
            "root_node": "left_subtree",
            "leaves": ["L1", "L2"],
            "leaf_signature": ("L1", "L2"),
            "size": 2,
        },
        1: {
            "root_node": "right_leaf",
            "leaves": ["R"],
            "leaf_signature": ("R",),
            "size": 1,
        },
    }


def test_build_sample_cluster_assignments_carries_leaf_signature() -> None:
    decomposition = {
        "cluster_assignments": {
            7: {
                "root_node": "internal_42",
                "leaves": ["S3", "S1"],
                "leaf_signature": ("S1", "S3"),
                "size": 2,
            }
        }
    }

    assignments = build_sample_cluster_assignments(decomposition)

    assert assignments.loc["S1", "cluster_root"] == "internal_42"
    assert assignments.loc["S1", "cluster_leaf_signature"] == ("S1", "S3")
    assert assignments.loc["S3", "cluster_leaf_signature"] == ("S1", "S3")


def test_build_sample_cluster_assignments_requires_cluster_assignments_key() -> None:
    with pytest.raises(KeyError, match="cluster_assignments"):
        build_sample_cluster_assignments({})


def test_build_sample_cluster_assignments_rejects_malformed_cluster_metadata() -> None:
    decomposition = {"cluster_assignments": {1: {"root_node": "left", "leaves": ["S1"]}}}

    with pytest.raises(KeyError, match="size"):
        build_sample_cluster_assignments(decomposition)


def test_build_sample_cluster_assignments_rejects_mismatched_leaf_signature() -> None:
    decomposition = {
        "cluster_assignments": {
            1: {
                "root_node": "left",
                "leaves": ["S1", "S3"],
                "leaf_signature": ("S1",),
                "size": 2,
            }
        }
    }

    with pytest.raises(ValueError, match="leaf_signature"):
        build_sample_cluster_assignments(decomposition)


def test_build_sample_cluster_assignments_rejects_overlapping_leaves() -> None:
    decomposition = {
        "cluster_assignments": {
            1: {"root_node": "left", "leaves": ["S1"], "size": 1},
            2: {"root_node": "right", "leaves": ["S1"], "size": 1},
        }
    }

    with pytest.raises(ValueError, match="multiple clusters"):
        build_sample_cluster_assignments(decomposition)
