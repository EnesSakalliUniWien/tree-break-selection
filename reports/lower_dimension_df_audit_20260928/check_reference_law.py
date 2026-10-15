"""Bounded Gaussian controls: dimension bookkeeping versus the null law."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

from benchmarks.diagnostics.calibration.statistics.statistic_distribution_shape_panel import (
    infer_satterthwaite_reference_from_eigenvalues,
)
from benchmarks.validation.contracts.report_contract import wilson_interval
from tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald.projected_wald_kernel import (
    run_projected_wald_kernel,
)

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
SEED, N, WIDTH, ALPHA = 2026092803, 100_000, 8, 0.01
SOURCES = [
    "tree_break_selection/eigenvalue_groups.py",
    "tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_kernel.py",
    "tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_projection_basis.py",
    "tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_reference_distribution.py",
    "benchmarks/diagnostics/calibration/statistics/statistic_distribution_shape_panel.py",
]


def hashes():
    return {path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest() for path in SOURCES}


def main():
    before = hashes()
    rng = np.random.default_rng(SEED)
    z = rng.normal(size=(N, WIDTH))
    independent = rng.normal(size=(N, WIDTH))
    independent /= np.linalg.norm(independent, axis=1, keepdims=True)
    energy = z * z
    total = energy.sum(axis=1)
    cumulative = np.cumsum(energy, axis=1)
    adaptive_k = np.argmax(cumulative >= 0.8 * total[:, None], axis=1) + 1
    adaptive_t = cumulative[np.arange(N), adaptive_k - 1]
    correlated = np.column_stack((z[:, 0], 0.8 * z[:, 0] + 0.6 * z[:, 1]))
    correlated_t = np.sum(correlated**2, axis=1)
    scale, effective_df = infer_satterthwaite_reference_from_eigenvalues([1.8, 0.2])
    cases = {
        "fixed_coordinate_rank1": (energy[:, 0], np.ones(N), 1.0),
        "independent_direction_rank1": ((independent * z).sum(axis=1) ** 2, np.ones(N), 1.0),
        "full_fixed_rank8": (total, np.full(N, WIDTH), 1.0),
        "contrast_aligned_rank1_nominal_df1": (total, np.ones(N), 1.0),
        "contrast_aligned_rank1_exact_df8": (total, np.full(N, WIDTH), 1.0),
        "contrast_energy_0.8_actual_selected_df": (adaptive_t, adaptive_k, 1.0),
        "correlated_fixed_rank2_nominal_df2": (correlated_t, np.full(N, 2), 1.0),
        "correlated_fixed_rank2_satterthwaite": (correlated_t, np.full(N, effective_df), scale),
        "correlated_fixed_rank2_known_covariance_whitened": (
            energy[:, :2].sum(axis=1),
            np.full(N, 2),
            1.0,
        ),
    }
    checked = 0
    for i in range(256):
        specifications = [
            (np.eye(WIDTH), np.arange(WIDTH, 0, -1), 1, None, energy[i, 0], 1),
            (
                independent[i : i + 1],
                np.ones(1),
                1,
                None,
                cases["independent_direction_rank1"][0][i],
                1,
            ),
            (z[i : i + 1] / np.sqrt(total[i]), np.ones(1), 1, None, total[i], 1),
            (np.eye(WIDTH), np.arange(WIDTH, 0, -1), WIDTH, 0.8, adaptive_t[i], adaptive_k[i]),
        ]
        for basis, eigenvalues, rank, fraction, expected_t, expected_df in specifications:
            result = run_projected_wald_kernel(
                z[i],
                spectral_k=rank,
                pca_projection=basis,
                pca_eigenvalues=eigenvalues,
                adaptive_dimension_energy_fraction=fraction,
            )
            np.testing.assert_allclose(result.statistic, expected_t, rtol=1e-12, atol=1e-12)
            assert result.degrees_of_freedom == expected_df == result.projection_dimension
            np.testing.assert_allclose(result.p_value, chi2.sf(expected_t, expected_df))
            checked += 1
    tie = run_projected_wald_kernel(
        z[0],
        spectral_k=1,
        pca_projection=np.eye(WIDTH),
        pca_eigenvalues=np.array([3, 3, 2, 1.8, 1.6, 1.4, 1.2, 1.0]),
    )
    assert tie.projection_dimension == tie.degrees_of_freedom == 2
    np.testing.assert_allclose(tie.statistic, energy[0, :2].sum())
    summaries = []
    for name, (statistic, df, reference_scale) in cases.items():
        p = chi2.sf(statistic / reference_scale, df)
        count = int((p <= ALPHA).sum())
        low, high = wilson_interval(count, N)
        summaries.append(
            {
                "case": name,
                "n": N,
                "alpha": ALPHA,
                "rejections": count,
                "rate": count / N,
                "wilson95_low": low,
                "wilson95_high": high,
                "df_min": float(df.min()),
                "df_max": float(df.max()),
                "reference_scale": reference_scale,
                "mean_statistic": statistic.mean(),
            }
        )
    pd.DataFrame(summaries).to_csv(BASE / "reference_law_summary.csv", index=False)
    after = hashes()
    assert before == after
    result = {
        "passed_arithmetic_checks": True,
        "seed": SEED,
        "independent_contrasts": N,
        "ambient_dimension": WIDTH,
        "alpha": ALPHA,
        "production_kernel_cross_checks": checked,
        "tie_extended_actual_df": 2,
        "aligned_rank1_exact_null_law": "chi_square_8 despite rank one",
        "aligned_rank1_exact_rejection_probability_using_df1": float(
            chi2.sf(chi2.isf(ALPHA, 1), WIDTH)
        ),
        "correlated_example_operator_weights": [1.8, 0.2],
        "correlated_example_satterthwaite_scale": scale,
        "correlated_example_satterthwaite_df": effective_df,
        "source_sha256_before_and_after": before,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope": "Independent Gaussian contrasts with known covariance, no tree or empirical inflation. Aligned-direction case is an analytic counterexample, not the production PCA algorithm. Energy-fraction case implements the production kernel rank rule with a fixed supplied basis. Satterthwaite is a moment approximation, not exact calibration. Method rows share matched draws.",
    }
    (BASE / "reference_law_verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
