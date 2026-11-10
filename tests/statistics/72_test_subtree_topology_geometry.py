import numpy as np
import pytest
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist

from benchmarks.diagnostics.calibration.selected.hierarchy.subtree_topology_geometry import (
    Cut,
    cophenetic_distances,
    fit_additive_distances,
    rooted_binary_topologies,
    rooted_robinson_foulds,
    spectral_geometry,
    spectral_subspaces,
    tree_clades,
    tree_edge_splits,
    unordered_cuts,
    ward_topology_family,
)


def test_cut_and_rooted_tree_families_have_different_counts():
    cuts = unordered_cuts(4)
    assert len(cuts) == len(set(cuts)) == 7
    assert all(0 in cut.left and set(cut.left).isdisjoint(cut.right) for cut in cuts)
    assert all(sorted(cut.left + cut.right) == list(range(4)) for cut in cuts)
    trees = rooted_binary_topologies(4)
    assert len(trees) == len(set(trees)) == 15
    assert len(rooted_binary_topologies(6)) == 945
    assert {cut for tree in trees for cut in tree_edge_splits(tree)} == set(cuts)
    for tree in trees:
        a, b = tree_clades(tree)
        assert a <= b or b <= a or not (a & b)
    left = ((0, 1), (2, 3))
    right = ((0, 2), (1, 3))
    assert rooted_robinson_foulds(left, left) == 0
    assert rooted_robinson_foulds(left, right) == 4
    assert len(tree_edge_splits(left)) == 5
    np.testing.assert_array_equal(cophenetic_distances(left), [1, 2, 2, 2, 2, 1])
    np.testing.assert_array_equal(cophenetic_distances(right), [2, 1, 2, 2, 1, 2])


def test_whole_tied_block_rotation_preserves_geometry_and_ward_tie_convention():
    square = np.array([[-1, -1], [-1, 1], [1, -1], [1, 1]], dtype=float)
    angle = 0.42
    rot = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    (first,) = spectral_geometry(square)
    (rotated,) = spectral_geometry(square @ rot)
    assert first.subspace.eigenblocks == rotated.subspace.eigenblocks == ((0, 2),)
    assert first.subspace.rank == rotated.subspace.rank == 2
    np.testing.assert_allclose(first.subspace.projector, rotated.subspace.projector, atol=1e-12)
    np.testing.assert_allclose(first.distances, rotated.distances)
    np.testing.assert_allclose(first.distances, pdist(square))
    assert first.cuts == rotated.cuts
    assert len(first.cuts) == 6 and Cut((0, 3), (1, 2)) not in first.cuts
    assert first.ward_tree == rotated.ward_tree
    np.testing.assert_allclose(
        cophenetic_distances(first.ward_linkage), cophenet(first.ward_linkage)
    )


def test_standalone_and_cumulative_blocks_and_saturated_span():
    tall = np.array([[-2, -1], [-2, 1], [2, -1], [2, 1]], dtype=float)
    spaces = {record.subspace.eigenblocks: record for record in spectral_geometry(tall)}
    assert {ranges: item.subspace.rank for ranges, item in spaces.items()} == {
        ((0, 1),): 1,
        ((1, 2),): 1,
        ((0, 1), (1, 2)): 2,
    }
    assert spaces[((0, 1),)].cuts == (Cut((0, 1), (2, 3)),)
    assert spaces[((1, 2),)].cuts == (Cut((0, 2), (1, 3)),)
    (saturated,) = spectral_geometry(np.eye(4))
    assert saturated.subspace.rank == 3 and saturated.cuts == unordered_cuts(4)
    np.testing.assert_allclose(saturated.subspace.projector, np.eye(4) - 0.25, atol=1e-12)


def test_coincident_positions_and_bounds():
    data = np.array([[0.0], [0.0], [2.0], [2.0]])
    (record,) = spectral_geometry(data)
    assert record.cuts == (Cut((0, 1), (2, 3)),)
    assert spectral_subspaces(np.ones((4, 2))) == ()
    assert spectral_subspaces(np.tile([0.1, 0.3, 1.0], (7, 1))) == ()
    assert spectral_geometry(np.ones((4, 2))) == ()
    with pytest.raises(ValueError, match="10"):
        spectral_geometry(np.ones((11, 2)))
    with pytest.raises(ValueError, match="6"):
        rooted_binary_topologies(7)
    with pytest.raises(ValueError, match="finite"):
        spectral_subspaces([[0.0], [np.inf]])


@pytest.mark.parametrize("scale", [1e-200, 1e200])
def test_distance_and_cut_calculation_preserve_uniform_feature_rescaling(scale):
    points = np.array([[0.0], [1.0], [3.0], [4.0]])
    (ordinary,) = spectral_geometry(points)
    (rescaled,) = spectral_geometry(points * scale)
    assert ordinary.cuts == rescaled.cuts
    assert ordinary.ward_tree == rescaled.ward_tree
    np.testing.assert_allclose(ordinary.distances, rescaled.distances / scale)


def test_ward_family_enumerates_and_deduplicates_all_tied_merges():
    vertical = ((0, 1), (2, 3))
    horizontal = ((0, 2), (1, 3))
    square = np.array([[-1, -1], [-1, 1], [1, -1], [1, 1]])
    assert ward_topology_family(square) == {vertical, horizontal}
    assert ward_topology_family(pdist(square)) == {vertical, horizontal}
    assert ward_topology_family(square[:, :1]) == {vertical}
    assert ward_topology_family(square[:, 1:]) == {horizontal}
    assert ward_topology_family(np.eye(4), rtol=0) == set(rooted_binary_topologies(4))
    assert ward_topology_family(np.eye(6)) == set(rooted_binary_topologies(6))
    (simplex,) = spectral_geometry(np.eye(4))
    assert simplex.ward_tree in ward_topology_family(simplex.distances)


def test_ward_family_tolerance_and_input_contract():
    almost_square = np.array([[-1.0, -1], [-1, 1], [1 + 1e-11, -1], [1 + 1e-11, 1]])
    assert len(ward_topology_family(almost_square, rtol=0)) == 1
    assert len(ward_topology_family(almost_square)) == 2
    assert len(ward_topology_family(np.array([[0], [1], [4], [10]]), rtol=0)) == 1
    with pytest.raises(ValueError, match="6"):
        ward_topology_family(np.zeros((7, 2)))
    with pytest.raises(ValueError, match="Euclidean"):
        ward_topology_family(np.array([1.0, 1.0, 3.0]))
    with pytest.raises(ValueError, match="rtol"):
        ward_topology_family(almost_square, rtol=-1)


def test_additive_distance_fit_recovers_lengths_and_is_invariant_to_root():
    observed = np.array([3, 8, 9, 9, 10, 7], dtype=float)
    balanced = fit_additive_distances(((0, 1), (2, 3)), observed)
    rerooted = fit_additive_distances((0, (1, (2, 3))), observed)
    wrong = fit_additive_distances(((0, 2), (1, 3)), observed)
    assert balanced.edge_splits == rerooted.edge_splits
    expected = {
        Cut((0,), (1, 2, 3)): 1,
        Cut((0, 2, 3), (1,)): 2,
        Cut((0, 1, 3), (2,)): 3,
        Cut((0, 1, 2), (3,)): 4,
        Cut((0, 1), (2, 3)): 4,
    }
    np.testing.assert_allclose(balanced.branch_lengths, [expected[s] for s in balanced.edge_splits])
    np.testing.assert_allclose(balanced.branch_lengths, rerooted.branch_lengths)
    np.testing.assert_allclose(balanced.fitted_distances, observed)
    assert balanced.rmse < 1e-12 and balanced.max_abs_residual < 1e-12
    assert wrong.rmse > 0.1 and wrong.max_abs_residual > 0.1
    assert wrong.branch_lengths[wrong.edge_splits.index(Cut((0, 2), (1, 3)))] == 0


def test_additive_distance_fit_input_contract():
    fit = fit_additive_distances((0, 1), np.array([0.0]))
    np.testing.assert_array_equal(fit.fitted_distances, [0])
    assert fit.rmse == fit.max_abs_residual == 0
    for tree in (((0, 0), 1), (0, 2)):
        with pytest.raises(ValueError, match="labels"):
            fit_additive_distances(tree, np.array([1.0]))
    for distances in ([1.0], [1.0, -2.0, 3.0], [1.0, np.inf, 3.0]):
        with pytest.raises(ValueError, match="distances"):
            fit_additive_distances(((0, 1), 2), distances)
