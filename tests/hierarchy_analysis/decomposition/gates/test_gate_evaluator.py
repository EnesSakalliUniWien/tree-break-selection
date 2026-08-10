"""Contracts for gate traversal decisions."""

from __future__ import annotations

import networkx as nx
import pytest
from tree_break_selection.hierarchy_analysis.decomposition.gates.gate_evaluator import (
    GateEvaluator,
    TraversalDecision,
)

from ...gate_support import (
    _annotate_tree_structure,
    _make_binary_tree,
    _make_deep_tree,
)


def _make_gate(
    tree: nx.DiGraph | None = None,
    *,
    edge_divergent: dict[str, bool] | None = None,
    sibling_different: dict[str, bool] | None = None,
    sibling_skipped: dict[str, bool] | None = None,
    children_map: dict[str, list[str]] | None = None,
    passthrough: bool = False,
    passthrough_supported: dict[str, bool] | None = None,
) -> GateEvaluator:
    if tree is None:
        tree = _make_binary_tree()

    if children_map is None:
        children_map = {node: list(tree.successors(node)) for node in tree.nodes}

    if edge_divergent is None:
        edge_divergent = {node: True for node in tree.nodes}

    if sibling_different is None:
        sibling_different = {node: True for node in tree.nodes}

    if sibling_skipped is None:
        sibling_skipped = {node: False for node in tree.nodes}

    return GateEvaluator(
        tree=tree,
        edge_divergent=edge_divergent,
        sibling_different=sibling_different,
        sibling_skipped=sibling_skipped,
        children_map=children_map,
        passthrough=passthrough,
        passthrough_supported=passthrough_supported,
    )


class TestGateEvaluator:
    def test_gate1_nonbinary_node_returns_false(self) -> None:
        tree = nx.DiGraph()
        tree.add_edges_from([("root", "A"), ("root", "B"), ("root", "C")])
        _annotate_tree_structure(tree, {"A", "B", "C"})

        gate = _make_gate(
            tree=tree,
            children_map={"root": ["A", "B", "C"], "A": [], "B": [], "C": []},
        )
        assert gate.decision("root") is TraversalDecision.BOUNDARY
        assert (
            gate.passthrough_audit_status("root")["passthrough_split_prerequisites_open"] is False
        )

    def test_gate1_single_child_returns_false(self) -> None:
        tree = nx.DiGraph()
        tree.add_edges_from([("root", "A")])
        _annotate_tree_structure(tree, {"A"})

        gate = _make_gate(
            tree=tree,
            children_map={"root": ["A"], "A": []},
        )
        assert gate.decision("root") is TraversalDecision.BOUNDARY
        assert (
            gate.passthrough_audit_status("root")["passthrough_split_prerequisites_open"] is False
        )

    def test_edge_gate_neither_child_diverges(self) -> None:
        gate = _make_gate(
            edge_divergent={
                "root": True,
                "L": False,
                "R": False,
                "L1": False,
                "L2": False,
                "R1": False,
                "R2": False,
            },
        )
        assert gate.decision("root") is TraversalDecision.BOUNDARY
        assert (
            gate.passthrough_audit_status("root")["passthrough_split_prerequisites_open"] is False
        )

    def test_edge_gate_one_child_diverges(self) -> None:
        gate = _make_gate(
            edge_divergent={
                "root": True,
                "L": True,
                "R": False,
                "L1": True,
                "L2": True,
                "R1": True,
                "R2": True,
            },
        )
        assert gate.decision("root") is TraversalDecision.SPLIT

    def test_edge_gate_missing_annotations_raises(self) -> None:
        with pytest.raises(ValueError, match="Missing edge_divergent values"):
            _make_gate(edge_divergent={})

    def test_sibling_gate_siblings_same(self) -> None:
        gate = _make_gate(
            sibling_different={
                "root": False,
                "L": False,
                "R": False,
                "L1": False,
                "L2": False,
                "R1": False,
                "R2": False,
            },
        )
        assert gate.decision("root") is TraversalDecision.BOUNDARY

    def test_sibling_gate_siblings_different(self) -> None:
        assert _make_gate().decision("root") is TraversalDecision.SPLIT

    def test_sibling_gate_skipped_returns_false(self) -> None:
        gate = _make_gate(
            sibling_skipped={
                "root": True,
                "L": False,
                "R": False,
                "L1": False,
                "L2": False,
                "R1": False,
                "R2": False,
            },
        )
        assert gate.decision("root") is TraversalDecision.BOUNDARY

    def test_sibling_gate_missing_annotations_raises(self) -> None:
        with pytest.raises(ValueError, match="Missing sibling_different values"):
            _make_gate(sibling_different={})

    def test_passthrough_disabled_returns_false(self) -> None:
        gate = _make_gate(
            passthrough=False,
            sibling_different={
                "root": False,
                "L": False,
                "R": False,
                "L1": False,
                "L2": False,
                "R1": False,
                "R2": False,
            },
        )
        assert gate.decision("root") is TraversalDecision.BOUNDARY
        assert gate.passthrough_audit_status("root") == {
            "passthrough_enabled": False,
            "passthrough_split_prerequisites_open": True,
            "passthrough_sibling_gate_open": False,
            "passthrough_descendant_split_available": False,
            "passthrough_candidate": False,
            "passthrough_supported": True,
            "passthrough_bottleneck": "",
            "passthrough_decision_reason": "passthrough_disabled",
        }

    def test_passthrough_when_sibling_gate_fails_with_descendant_signal(self) -> None:
        tree = _make_deep_tree()
        sibling_different = {node: False for node in tree.nodes}
        sibling_different["B"] = True
        gate = _make_gate(
            tree=tree,
            passthrough=True,
            sibling_different=sibling_different,
        )
        assert gate.decision("root") is TraversalDecision.PASS_THROUGH
        assert gate.passthrough_audit_status("root") == {
            "passthrough_enabled": True,
            "passthrough_split_prerequisites_open": True,
            "passthrough_sibling_gate_open": False,
            "passthrough_descendant_split_available": True,
            "passthrough_candidate": True,
            "passthrough_supported": True,
            "passthrough_bottleneck": "",
            "passthrough_decision_reason": "pass_through",
        }

    def test_spectral_support_guard_blocks_passthrough(self) -> None:
        tree = _make_deep_tree()
        sibling_different = {node: False for node in tree.nodes}
        sibling_different["B"] = True
        passthrough_supported = {node: True for node in tree.nodes}
        passthrough_supported["root"] = False

        gate = _make_gate(
            tree=tree,
            passthrough=True,
            sibling_different=sibling_different,
            passthrough_supported=passthrough_supported,
        )

        assert gate.decision("root") is TraversalDecision.BOUNDARY
        assert gate.passthrough_support_status("root") == {
            "passthrough_supported": False,
            "passthrough_bottleneck": "",
        }
        assert gate.passthrough_audit_status("root") == {
            "passthrough_enabled": True,
            "passthrough_split_prerequisites_open": True,
            "passthrough_sibling_gate_open": False,
            "passthrough_descendant_split_available": True,
            "passthrough_candidate": True,
            "passthrough_supported": False,
            "passthrough_bottleneck": "",
            "passthrough_decision_reason": "passthrough_support_blocked",
        }

    def test_passthrough_decision_uses_cached_gate_results(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        tree = _make_deep_tree()
        sibling_different = {node: False for node in tree.nodes}
        sibling_different["B"] = True
        gate = _make_gate(
            tree=tree,
            passthrough=True,
            sibling_different=sibling_different,
        )

        def _fail_recompute(*_args, **_kwargs):
            raise AssertionError("passthrough decision recomputed gate predicates")

        monkeypatch.setattr(gate, "_passes_split_prerequisites", _fail_recompute)
        monkeypatch.setattr(gate, "_sibling_gate_is_open", _fail_recompute)

        assert gate.decision("root") is TraversalDecision.PASS_THROUGH

    def test_no_passthrough_when_sibling_gate_passes(self) -> None:
        gate = _make_gate(
            passthrough=True,
        )
        assert gate.decision("root") is TraversalDecision.SPLIT
        assert (
            gate.passthrough_audit_status("root")["passthrough_decision_reason"]
            == "sibling_gate_open_split"
        )

    def test_no_passthrough_when_gates_1_2_fail(self) -> None:
        gate = _make_gate(
            passthrough=True,
            edge_divergent={
                "root": False,
                "L": False,
                "R": False,
                "L1": False,
                "L2": False,
                "R1": False,
                "R2": False,
            },
        )
        assert gate.decision("root") is TraversalDecision.BOUNDARY
        assert (
            gate.passthrough_audit_status("root")["passthrough_decision_reason"]
            == "split_prerequisites_closed"
        )

    def test_no_passthrough_when_no_descendant_signal(self) -> None:
        gate = _make_gate(
            passthrough=True,
            sibling_different={
                "root": False,
                "L": False,
                "R": False,
                "L1": False,
                "L2": False,
                "R1": False,
                "R2": False,
            },
        )
        assert gate.decision("root") is TraversalDecision.BOUNDARY
        assert (
            gate.passthrough_audit_status("root")["passthrough_decision_reason"]
            == "no_descendant_split"
        )
