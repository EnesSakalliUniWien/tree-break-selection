from __future__ import annotations

import pytest
from benchmarks.diagnostics.calibration.selected.hierarchy.selected_hierarchy_null_audit import (
    SelectedSiblingRecord,
    _choose_observed_target,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.pair_testing.types.sibling_pair_record import (
    SiblingPairRecord,
)


def _selected(parent: str, *, depth: int, n_parent: int, p_value: float) -> SelectedSiblingRecord:
    record = SiblingPairRecord(
        parent=parent,
        left=f"{parent}L",
        right=f"{parent}R",
        stat=4.0,
        reference_scale=1.0,
        degrees_of_freedom=2.0,
        p_value=p_value,
        branch_length_sum=0.1,
        n_parent=n_parent,
        is_null_like=False,
        sibling_projection_dimension=2.0,
        feature_family="bernoulli",
    )
    return SelectedSiblingRecord(record=record, parent_depth=depth)


def test_root_target_is_the_depth_zero_parent() -> None:
    records = (
        _selected("root", depth=0, n_parent=100, p_value=0.9),
        _selected("N1", depth=1, n_parent=60, p_value=0.1),
    )
    assert _choose_observed_target("case", records, target_mode="root").parent == "root"
    assert _choose_observed_target("case", records, target_mode="non_root_strongest").parent == "N1"


def test_root_target_is_not_substituted_by_the_largest_open_parent() -> None:
    # The root edge was blocked, so only non-root records survive the edge path.
    # N1 is the largest open parent; the old proxy mislabelled it root and dropped it.
    records = (
        _selected("N1", depth=1, n_parent=60, p_value=0.1),
        _selected("N2", depth=2, n_parent=30, p_value=0.5),
    )
    with pytest.raises(ValueError, match="no root sibling record"):
        _choose_observed_target("case", records, target_mode="root")
    assert _choose_observed_target("case", records, target_mode="non_root_strongest").parent == "N1"
