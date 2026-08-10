"""Contracts for spectral-transport passthrough annotation."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
from tree_break_selection.hierarchy_analysis.decomposition.gates.spectral_transport import (
    annotate_spectral_transport_passthrough_support,
)

from ...gate_support import (
    _make_annotations,
    _make_deep_tree,
)


def test_spectral_transport_annotation_supports_coherent_passthrough_path() -> None:
    tree = _make_deep_tree()
    edge_divergent = {node: True for node in tree.nodes}
    sibling_different = {node: False for node in tree.nodes}
    sibling_different["B"] = True
    annotations = _make_annotations(
        tree,
        edge_divergent=edge_divergent,
        sibling_different=sibling_different,
    )
    spectral_context = SimpleNamespace(
        principal_component_projections_by_node={
            node: np.asarray([[1.0, 0.0]]) for node in tree.nodes
        },
        principal_component_eigenvalues_by_node={node: np.asarray([4.0]) for node in tree.nodes},
        raw_mp_signal_counts_by_node={node: 1 for node in tree.nodes},
    )

    out = annotate_spectral_transport_passthrough_support(
        tree,
        annotations,
        spectral_context,
        max_cost=0.1,
    )

    assert bool(out.loc["root", "Spectral_Transport_Pass_Through_Supported"])
    assert out.loc["root", "Spectral_Transport_Bottleneck"] == "supported_mp_mode_path"


def test_spectral_transport_annotation_blocks_rotated_passthrough_path() -> None:
    tree = _make_deep_tree()
    edge_divergent = {node: True for node in tree.nodes}
    sibling_different = {node: False for node in tree.nodes}
    sibling_different["B"] = True
    annotations = _make_annotations(
        tree,
        edge_divergent=edge_divergent,
        sibling_different=sibling_different,
    )
    spectral_context = SimpleNamespace(
        principal_component_projections_by_node={
            **{node: np.asarray([[1.0, 0.0]]) for node in tree.nodes},
            "B": np.asarray([[0.0, 1.0]]),
        },
        principal_component_eigenvalues_by_node={node: np.asarray([4.0]) for node in tree.nodes},
        raw_mp_signal_counts_by_node={node: 1 for node in tree.nodes},
    )

    out = annotate_spectral_transport_passthrough_support(
        tree,
        annotations,
        spectral_context,
        max_cost=0.1,
    )

    assert not bool(out.loc["root", "Spectral_Transport_Pass_Through_Supported"])
    assert bool(out.loc["root", "Spectral_Transport_Pass_Through_Blocked"])


def test_spectral_transport_annotation_blocks_floor_only_when_mp_blocks_required() -> None:
    tree = _make_deep_tree()
    edge_divergent = {node: True for node in tree.nodes}
    sibling_different = {node: False for node in tree.nodes}
    sibling_different["B"] = True
    annotations = _make_annotations(
        tree,
        edge_divergent=edge_divergent,
        sibling_different=sibling_different,
    )
    spectral_context = SimpleNamespace(
        principal_component_projections_by_node={
            node: np.asarray([[1.0, 0.0]]) for node in tree.nodes
        },
        principal_component_eigenvalues_by_node={node: np.asarray([4.0]) for node in tree.nodes},
        raw_mp_signal_counts_by_node={node: 0 for node in tree.nodes},
    )

    out = annotate_spectral_transport_passthrough_support(
        tree,
        annotations,
        spectral_context,
        max_cost=0.1,
    )

    assert not bool(out.loc["root", "Spectral_Transport_Pass_Through_Supported"])
    assert bool(out.loc["root", "Spectral_Transport_Pass_Through_Blocked"])
    assert out.loc["root", "Spectral_Transport_Bottleneck"] == ("spectral_transport_bottleneck")


def test_spectral_transport_annotation_does_not_veto_floor_only_when_mp_blocks_not_required() -> (
    None
):
    tree = _make_deep_tree()
    edge_divergent = {node: True for node in tree.nodes}
    sibling_different = {node: False for node in tree.nodes}
    sibling_different["B"] = True
    annotations = _make_annotations(
        tree,
        edge_divergent=edge_divergent,
        sibling_different=sibling_different,
    )
    spectral_context = SimpleNamespace(
        principal_component_projections_by_node={
            node: np.asarray([[1.0, 0.0]]) for node in tree.nodes
        },
        principal_component_eigenvalues_by_node={node: np.asarray([4.0]) for node in tree.nodes},
        raw_mp_signal_counts_by_node={node: 0 for node in tree.nodes},
    )

    out = annotate_spectral_transport_passthrough_support(
        tree,
        annotations,
        spectral_context,
        max_cost=0.1,
        require_mp_blocks=False,
    )

    assert bool(out.loc["root", "Spectral_Transport_Pass_Through_Supported"])
    assert not bool(out.loc["root", "Spectral_Transport_Pass_Through_Blocked"])
    assert out.loc["root", "Spectral_Transport_Bottleneck"] == ("unmeasured_no_matched_mp_path")
