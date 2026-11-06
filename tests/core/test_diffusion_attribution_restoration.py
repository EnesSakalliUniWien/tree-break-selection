"""Crossed attribution preserves the fixed tree and declared edge-role population."""

import numpy as np
import pytest

from benchmarks.diagnostics.calibration.diffusion_attribution import (
    full_space_sibling_records,
    gaussian_tree_pairs,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_attribution_restoration import (
    ROLE_FIELDS,
    compare_matrix,
    write_panel,
)


def test_crossed_restoration_preserves_baseline_and_roles(tmp_path):
    pair = next(gaussian_tree_pairs(replicates=1))
    rows = compare_matrix(pair.tree, pair.independent_data, pair.feature_space)
    assert len(rows) == 79 * 2 * 2 * 2
    assert rows.groupby("parent").size().eq(8).all()
    assert rows.branch_length_sum.eq(0).all()
    projected = rows[rows.projection == "production"]
    assert projected.groupby("parent").nominal_projection_dimension.nunique().eq(1).all()
    assert projected.spectral_covariance_chart.eq("native_parent_total").all()
    withheld = rows[(rows.inflation == "empirical") & ~rows.gate_eligible]
    assert withheld.production_gate_p.isna().all()
    assert withheld.gate_status.eq("skipped_null_like").all()
    assert rows.edge_roles_recomputed.eq(False).all()
    for field in ROLE_FIELDS:
        assert rows.groupby("parent")[field].nunique(dropna=False).eq(1).all()
    _, baseline = full_space_sibling_records(
        pair.tree, pair.independent_data, pair.feature_space, "known_identity"
    )
    direct = rows[
        (rows.covariance == "known_identity")
        & (rows.projection == "full_feature")
        & (rows.inflation == "none")
    ].set_index("parent")
    for record in baseline:
        assert direct.loc[record.parent, "raw_statistic"] == record.stat
        assert direct.loc[record.parent, "raw_df"] == record.degrees_of_freedom
        assert direct.loc[record.parent, "adjusted_p"] == record.p_value
    for _, group in rows.groupby(["parent", "covariance", "projection"]):
        assert group.raw_statistic.nunique() == 1
        assert group.raw_df.nunique() == 1
        empirical = group[group.inflation == "empirical"].iloc[0]
        if empirical.available:
            assert empirical.c_hat >= 1.0
            assert empirical.adjusted_p >= empirical.raw_p - 1e-15
        else:
            assert np.isnan(empirical.adjusted_p)
            assert empirical.calibration_status.startswith("undefined_")
    full_span = rows[(rows.projection == "production") & rows.full_span]
    for _, row in full_span.iterrows():
        reference = rows[
            (rows.parent == row.parent)
            & (rows.covariance == row.covariance)
            & (rows.projection == "full_feature")
            & (rows.inflation == row.inflation)
        ].iloc[0]
        assert row.raw_statistic == pytest.approx(reference.raw_statistic, rel=1e-10, abs=1e-10)
    rows["replicate"] = 0
    rows["data_source"] = "independent"
    summary = write_panel(rows, tmp_path)
    assert summary.independent_dataset_count.eq(1).all()
    assert (tmp_path / "restoration_contract.json").is_file()
    assert (tmp_path / "restoration_roots.csv").is_file()
