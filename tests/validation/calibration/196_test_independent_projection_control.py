import json
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import chi2

from benchmarks.diagnostics.calibration.sibling.independent_projection_control import (
    ControlConfig,
    draw_data,
    evaluate_projection,
    freeze,
    paired_outcomes,
    run_control,
    summarize,
    verify_frozen,
)
from benchmarks.diagnostics.calibration.sibling.known_branch_skew_diagnostic import (
    KnownBranchConfig,
    covariance_matrix,
    population_means,
    run_diagnostic,
)
from tree_break_selection.hierarchy_analysis.decomposition.core.eigen_result import EigenResult


def test_exact_evaluation_matching_and_separate_training():
    cfg = ControlConfig(
        replicates=2,
        n_per_group=8,
        dimension=4,
        patterns=("global_null",),
        covariance_profiles=("identity",),
    )
    rows, datasets = run_control(cfg)
    original = run_diagnostic(
        KnownBranchConfig(
            output_dir=Path("unused"),
            replicates=2,
            n_per_group=8,
            dimension=4,
            patterns=cfg.patterns,
            covariance_profiles=cfg.covariance_profiles,
        )
    )
    assert (
        datasets.dataset_sha256.tolist()
        == original["local"].dataset_sha256.drop_duplicates().tolist()
    )
    assert len(rows) == 2 * 3 * (15 + 6)
    for row in datasets.itertuples():
        assert (
            len({row.evaluation_sha256, row.basis_training_sha256, row.rank_training_sha256}) == 3
        )
        assert json.loads(row.basis_training_spawn_key)[1] == 1
        assert json.loads(row.rank_training_spawn_key)[1] == 2
    a = rows[
        (rows.regime == "canonical_mp")
        & (rows.basis_source == "evaluation")
        & (rows.rank_source == "evaluation_mp")
    ]
    b = original["local"].query(
        "covariance_estimator == 'oracle_covariance' and projection_kind == 'same_data_pca_mp'"
    )
    np.testing.assert_allclose(a.p_value, b.p_value, atol=1e-13)
    assert rows.status.eq("ok").all()


def test_actual_rank_extends_ties_but_fixed_coordinates_do_not():
    eigen = EigenResult(
        eigenvalues=np.array([3.0, 3.0, 1.0]),
        is_active_feature=np.ones(3, dtype=bool),
        active_feature_count=3,
        use_dual=False,
        eigenvectors_active=np.eye(3),
    )
    result = evaluate_projection(np.array([1.0, 2.0, 4.0]), eigen, 1)
    assert result["actual_rank"] == 2
    assert result["statistic"] == 5
    assert result["p_value"] == chi2.sf(5, 2)
    fixed = evaluate_projection(np.array([1.0, 2.0, 4.0]), None, 1)
    assert fixed["actual_rank"] == 1
    assert fixed["statistic"] == 1


def test_summaries_preserve_unavailable_and_paired_denominators():
    cfg = ControlConfig(
        replicates=2,
        n_per_group=8,
        dimension=4,
        patterns=("global_null",),
        covariance_profiles=("identity",),
    )
    rows, _ = run_control(cfg)
    rows.loc[0, ["status", "p_value"]] = ["runtime_error", np.nan]
    summary = summarize(rows, 0.01)
    assert summary.n_unavailable.sum() == 1
    paired, stats = paired_outcomes(rows, 0.01)
    assert paired.pair_status.eq("unavailable").any()
    assert stats.n_unavailable.sum() > 0
    assert (
        stats.both_reject + stats.a_only_reject + stats.b_only_reject + stats.neither_reject
    ).equals(stats.n_pairs)
    assert summary.query(
        "basis_source == 'independent' and rank_source == 'independent_mp'"
    ).conditional_chi_square_guarantee.all()


def test_sampling_is_repeatable_without_global_rng_state():
    means, sigma = population_means("global_null", 4), covariance_matrix("identity", 4)
    first = draw_data(np.random.SeedSequence(17, spawn_key=(1,)), means, sigma, 8)
    np.random.default_rng(45).standard_normal(100)
    second = draw_data(np.random.SeedSequence(17, spawn_key=(1,)), means, sigma, 8)
    np.testing.assert_array_equal(first, second)


def test_frozen_config_and_output_are_guarded(tmp_path):
    protocol = Path(
        "reports/independent_projection_vandermonde_20260928/control_protocol.md"
    ).resolve()
    cfg = ControlConfig(replicates=2)
    freeze(tmp_path, cfg, protocol, "test")
    verify_frozen(tmp_path, cfg, protocol)
    with pytest.raises(ValueError, match="Configuration"):
        verify_frozen(tmp_path, ControlConfig(replicates=3), protocol)
    (tmp_path / "result.csv").write_text("keep")
    with pytest.raises(FileExistsError):
        verify_frozen(tmp_path, cfg, protocol)
    with pytest.raises(FileExistsError):
        freeze(tmp_path, cfg, protocol, "test")
