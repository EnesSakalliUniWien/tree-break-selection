"""Descendant chain tests: ancestor-path events, node paths and the root-plus-children family."""

import numpy as np
import pytest

from benchmarks.diagnostics.calibration.selected.hierarchy import (
    diffusion_root_conditional_mc as mc,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import (
    make_diffusion_root_path,
)


def _hierarchical(rows_per_group, shift_top, shift_child, seed):
    values = np.random.default_rng(seed).normal(size=(4 * rows_per_group, 3))
    half = 2 * rows_per_group
    values[:half, 0] -= shift_top
    values[half:, 0] += shift_top
    for start in (0, half):
        values[start : start + rows_per_group, 1] -= shift_child
        values[start + rows_per_group : start + half, 1] += shift_child
    return values


def test_ancestor_path_event_requires_every_split_down_to_the_node():
    # Leaves 0..3; node 4 = (0, 1), node 5 = (2, 3), root 6 = (4, 5).
    observed = ((0, 1), (2, 3), (4, 5))
    required = mc._ancestor_path_pairs(observed, 4, 4)
    assert required == {((0,), (1,)), ((0, 1), (2, 3))}
    assert required <= mc._merge_pairs(((2, 3), (0, 1), (4, 5)), 4)
    assert not required <= mc._merge_pairs(((0, 2), (1, 3), (4, 5)), 4)
    # The root's own event is only its split.
    assert mc._ancestor_path_pairs(observed, 4, 6) == {((0, 1), (2, 3))}


def test_node_path_moves_only_the_node_rows_and_reproduces_the_data():
    path = make_diffusion_root_path(_hierarchical(4, 4.0, 2.0, seed=3), sigma=1.0)
    n = len(path.values)
    child = next(node for node in path.observed_signature[-1] if node >= n)
    node_path = mc._node_path(path, child)
    rows = mc._node_rows(path.observed_signature, n, child)
    outside = sorted(set(range(n)) - set(rows))

    np.testing.assert_allclose(
        node_path.offset + node_path.observed_radius * node_path.slope, path.values, atol=1e-12
    )
    np.testing.assert_array_equal(node_path.slope[outside], 0.0)
    np.testing.assert_allclose(node_path.e @ node_path.offset, 0.0, atol=1e-12)


def test_root_position_reproduces_the_validated_root_chain_test():
    values = np.random.default_rng(2026110401).normal(size=(9, 3))
    values[:4, 0] -= 2
    values[4:, 0] += 2
    root = mc.selected_diffusion_root_chain_test(values, sigma=1, steps=99, seed=7, alpha=0.05)
    family = mc.selected_diffusion_family_chain_tests(values, sigma=1, steps=99, seed=7, alpha=0.05)

    first = family["nodes"][0]
    assert first["role"] == "root"
    assert first["chain_radii"] == root["chain_radii"]
    assert first["p_value"] == root["p_value"]
    assert first["accepted_moves"] == root["accepted_moves"]


def test_family_is_the_root_and_its_internal_children_at_a_fixed_bonferroni_level():
    values = _hierarchical(3, 4.0, 2.0, seed=5)
    family = mc.selected_diffusion_family_chain_tests(
        values, sigma=1, steps=59, seed=11, alpha=0.05
    )

    assert family["selection_event"] == "ancestor_path"
    assert family["family_alpha"] == pytest.approx(0.05 / 3)
    roles = [node["role"] for node in family["nodes"]]
    assert roles[0] == "root" and set(roles[1:]) <= {"child"} and len(roles) <= 3
    for node in family["nodes"]:
        expected = (
            node["status"] == "complete"
            and node["sampling_status"] == "moved"
            and node["p_value"] <= family["family_alpha"]
        )
        assert (node["decision"] == "reject") == expected
        assert len(node["chain_radii"]) == 60


def test_strong_hierarchical_signal_rejects_the_root_and_both_children():
    # 32 rows: with fewer, the diffusion tree splits each child 2/6 across the subgroups.
    values = _hierarchical(8, 8.0, 6.0, seed=7)
    truth = np.repeat(np.arange(4), 8)
    family = mc.selected_diffusion_family_chain_tests(
        values, sigma=1, steps=999, seed=2026110600, alpha=0.05
    )

    assert [node["role"] for node in family["nodes"]] == ["root", "child", "child"]
    for node in family["nodes"][1:]:
        assert all(
            len({truth[int(label[1:])] for label in child}) == 1 for child in node["children"]
        )
    assert all(node["decision"] == "reject" for node in family["nodes"])


def test_family_reports_an_unsupported_observation_without_decisions():
    family = mc.selected_diffusion_family_chain_tests(np.zeros((9, 1)), sigma=1, steps=3)
    assert family["status"] == "unresolved" and family["nodes"] == []
    assert family["error_phase"] == "observation"
