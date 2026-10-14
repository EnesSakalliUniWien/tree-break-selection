"""Compare nominal, independently fitted, and exact Gaussian example laws."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from benchmarks.diagnostics.calibration.reporting import write_diagnostic_bundle
from benchmarks.validation.contracts.report_contract import wilson_interval
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.stats import chi2, gamma
from tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald.projected_wald_kernel import (
    run_projected_wald_kernel,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.selected_gaussian_radial import (
    RadialSelectionInterval,
    compute_selected_gaussian_tail,
)

matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
TRAIN_SEED, TEST_SEED = 2026092806, 2026092807
N_TRAIN, N_TEST = 100_000, 200_000
ALPHAS = (0.001, 0.005, 0.01, 0.05)
CASE_NAMES = (
    "fixed_rank1", "correlated_rank2", "energy_0.9_rank1",
    "energy_0.9_rank2", "energy_0.5_rank1", "radial_selected_rank2", "aligned_rank1",
)
LABELS = {
    "fixed_rank1": "Fixed rank one",
    "correlated_rank2": "Correlated fixed rank two",
    "energy_0.9_rank1": "Energy 0.9, selected rank one",
    "energy_0.9_rank2": "Energy 0.9, selected rank two",
    "energy_0.5_rank1": "Energy 0.5, selected rank one",
    "radial_selected_rank2": "Rank two, conditioned on T > 4",
    "aligned_rank1": "Eight-dimensional contrast, aligned rank one",
}
SOURCES = (
    "tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_kernel.py",
    "tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_projection_basis.py",
    "tree_break_selection/hierarchy_analysis/statistics/projection/projected_wald/projected_wald_reference_distribution.py",
    "tree_break_selection/hierarchy_analysis/statistics/projection/selected_gaussian_radial.py",
    "raw/inbox/wolfram-distribution-fitting-20260928.md",
    "raw/assets/benchmark-results/selected_tail_equation_cloud_run_20260604_1000/candidate_equation_holdout.csv",
    "raw/assets/benchmark-results/selected_tail_equation_cloud_run_20260604_1000/selected_ratio_tail_law.csv",
)


def hashes():
    return {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in SOURCES}


def generate(seed, n):
    z = np.random.default_rng(seed).standard_normal((n, 8))
    x, y = z[:, 0]**2, z[:, 1]**2
    s = x+y
    ratio = x/s
    return z, {
        "fixed_rank1": x,
        "correlated_rank2": 1.8*x + 0.2*y,
        "energy_0.9_rank1": x[ratio >= 0.9],
        "energy_0.9_rank2": s[ratio < 0.9],
        "energy_0.5_rank1": x[ratio >= 0.5],
        "radial_selected_rank2": s[s > 4],
        "aligned_rank1": (z*z).sum(axis=1),
    }


def angular_parameters(case):
    if case == "correlated_rank2":
        return np.pi/2, lambda angle: 1.8*np.cos(angle)**2 + 0.2*np.sin(angle)**2
    q = 0.9 if case == "energy_0.9_rank1" else 0.5
    return np.arccos(np.sqrt(q)), lambda angle: np.cos(angle)**2


def exact_sf(case, values, order=64):
    t = np.atleast_1d(np.asarray(values, dtype=float))
    if case in ("fixed_rank1", "energy_0.9_rank2", "aligned_rank1"):
        return chi2.sf(t, {"fixed_rank1": 1, "energy_0.9_rank2": 2, "aligned_rank1": 8}[case])
    if case == "radial_selected_rank2":
        return np.exp(-np.maximum(t-4, 0)/2)
    bound, weight = angular_parameters(case)
    nodes, weights = leggauss(order)
    factors = weight((nodes+1)*bound/2)
    output = np.empty(len(t))
    for start in range(0, len(t), 4096):
        chunk = t[start:start+4096]
        output[start:start+len(chunk)] = np.exp(-chunk[:, None]/(2*factors)) @ (weights/2)
    return output


def ks_uniform(p):
    ordered = np.sort(p)
    n = len(ordered)
    return max(float(np.max(np.arange(1, n+1)/n-ordered)), float(np.max(ordered-np.arange(n)/n)))


def main():
    source_hashes = hashes()
    design = {
        "train_seed": TRAIN_SEED, "test_seed": TEST_SEED,
        "n_train_unconditional": N_TRAIN, "n_test_unconditional": N_TEST,
        "alphas": ALPHAS, "cases": CASE_NAMES, "source_sha256": source_hashes,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "fit_policy": "Fit once to independent training contrasts, separately per case/selected rank. No tuning on evaluation data. Conditional samples have random sizes; cases share matched draws.",
    }
    (BASE/"pre_run_design.json").write_text(json.dumps(design, indent=2)+'\n')
    _, train = generate(TRAIN_SEED, N_TRAIN)
    z, test = generate(TEST_SEED, N_TEST)
    # Independent Wolfram Language numerical checkpoints, 35-digit integration.
    checkpoints = {
        "energy_0.9_rank1": [0.5957907655987167, 0.12620079953583704, 0.02679366463202124, 0.005701190670614923, 0.00025975594593373127],
        "correlated_rank2": [0.51685399553437, 0.14726102752233234, 0.052209612639442566, 0.01972154281568695, 0.003061821726528723],
    }
    max_quad_error = 0.0
    for case in ("energy_0.9_rank1", "energy_0.5_rank1", "correlated_rank2"):
        values = np.array([0., 1., 4., 7., 10., 16., 30., 60.])
        bound, weight = angular_parameters(case)
        reference = np.array([quad(lambda angle: np.exp(-t/(2*weight(angle)))/bound, 0, bound, epsabs=1e-12)[0] for t in values])
        error = float(np.max(np.abs(reference-exact_sf(case, values))))
        max_quad_error = max(max_quad_error, error)
        assert error < 1e-11
    for case, expected in checkpoints.items():
        np.testing.assert_allclose(exact_sf(case, [1,4,7,10,16]), expected, atol=1e-12, rtol=1e-11)
    for row in z[:256]:
        result = run_projected_wald_kernel(
            row[:2], spectral_k=2, pca_projection=np.eye(2), pca_eigenvalues=np.array([2.,1.]),
            adaptive_dimension_energy_fraction=0.9,
        )
        rank = 1 if row[0]**2 >= 0.9*np.sum(row[:2]**2) else 2
        assert result.projection_dimension == result.degrees_of_freedom == rank
        np.testing.assert_allclose(result.statistic, np.sum(row[:rank]**2))
    for t in (4.1, 7., 10., 16.):
        result = compute_selected_gaussian_tail(
            observed_radius=np.sqrt(t), projection_dimension=2,
            selection_intervals=(RadialSelectionInterval(2., np.inf),),
        )
        np.testing.assert_allclose(result.p_value, exact_sf("radial_selected_rank2", [t])[0], rtol=1e-12)
    parameter_rows, comparison_rows, curves = [], [], []
    plot_data = {}
    for case in CASE_NAMES:
        training, evaluation = train[case], test[case]
        mean, variance = float(training.mean()), float(training.var(ddof=1))
        current_df = 2 if case in ("correlated_rank2", "energy_0.9_rank2", "radial_selected_rank2") else 1
        shape, _, gamma_scale = gamma.fit(training, floc=0)
        parameters = {
            "current_chi_square": (1., float(current_df)),
            "df_only_moment_fit": (1., mean),
            "scale_only_moment_fit": (mean/current_df, float(current_df)),
            "scale_and_df_moment_fit": (variance/(2*mean), 2*mean**2/variance),
            "scale_and_df_mle_fit": (float(gamma_scale)/2, 2*float(shape)),
        }
        grid = np.linspace(0, max(15., float(np.quantile(evaluation, .999))), 301)
        true_curve = exact_sf(case, grid)
        all_p = {method: chi2.sf(evaluation/a, nu) for method, (a,nu) in parameters.items()}
        all_p["exact_conditional_law"] = exact_sf(case, evaluation)
        for method, p in all_p.items():
            if method in parameters:
                a, nu = parameters[method]
                parameter_rows.append({"case":case,"method":method,"n_train":len(training),"scale":a,"df":nu,"training_mean":mean,"training_variance":variance})
                curve = chi2.sf(grid/a, nu)
            else:
                curve = true_curve
            for t, truth, predicted in zip(grid, true_curve, curve, strict=True):
                curves.append({"case":case,"method":method,"statistic":t,"exact_survival":truth,"reference_survival":predicted})
            for alpha in ALPHAS:
                count = int((p <= alpha).sum())
                low, high = wilson_interval(count, len(evaluation))
                true_rate = float(exact_sf(case, [a*chi2.isf(alpha, nu)])[0]) if method in parameters else alpha
                comparison_rows.append({
                    "case":case,"method":method,"n_train":len(training),"n_test":len(evaluation),
                    "alpha":alpha,"rejections":count,"rate":count/len(evaluation),
                    "wilson95_low":low,"wilson95_high":high,"true_rate_for_frozen_fit":true_rate,
                    "true_rate_over_alpha":true_rate/alpha,"uniform_ks_distance":ks_uniform(p),
                    "cdf_grid_max_error":float(np.max(np.abs(curve-true_curve))),
                })
        plot_data[case] = (evaluation, grid, parameters)
    comparisons = pd.DataFrame(comparison_rows)
    mixture_rows = []
    rank1_probability = 2*np.arccos(np.sqrt(.9))/np.pi
    selected_rows = comparisons[comparisons.case.isin(["energy_0.9_rank1", "energy_0.9_rank2"])]
    for (method, alpha), group in selected_rows.groupby(["method", "alpha"]):
        n = int(group.n_test.sum())
        count = int(group.rejections.sum())
        assert n == N_TEST
        low, high = wilson_interval(count, n)
        true_rate = sum(
            (rank1_probability if row.case == "energy_0.9_rank1" else 1-rank1_probability)
            * row.true_rate_for_frozen_fit for row in group.itertuples()
        )
        mixture_rows.append({"method":method,"alpha":alpha,"n_test":n,"rejections":count,
                             "rate":count/n,"wilson95_low":low,"wilson95_high":high,
                             "true_rate_for_frozen_fit":true_rate})
    historical = pd.read_csv(ROOT/SOURCES[-2])
    historical = historical[historical.split_strategy.eq("leave_one_source_family_out")]
    tail_contexts = pd.read_csv(ROOT/SOURCES[-1])
    assert len(historical) == 6
    assert hashes() == source_hashes
    manifest = {
        **design, "scope":"Synthetic Gaussian examples, known covariance and declared selection events. No full tree inference, same-data learned PCA basis, fitted empirical sibling inflation, categorical test or production change.",
        "wolfram_numeric_checkpoints":10,"quadrature_max_absolute_error":max_quad_error,
        "actual_production_kernel_checks":256,"existing_selected_radial_tail_checks":4,
        "source_hashes_unchanged":True,
        "historical_tail_contexts":len(tail_contexts),
        "historically_marked_admissible_contexts":int(tail_contexts.production_tail_law_admissible.sum()),
    }
    write_diagnostic_bundle(
        output_dir=BASE, tables={"parameters":pd.DataFrame(parameter_rows),"comparison":comparisons,"curves":pd.DataFrame(curves),"rank_mixture":pd.DataFrame(mixture_rows),"historical_equations":historical},
        filenames={"parameters":"fitted_parameters.csv","comparison":"heldout_comparison.csv","curves":"distribution_curves.csv","rank_mixture":"rank_mixture_summary.csv","historical_equations":"historical_equation_holdout.csv"},manifest=manifest,
    )
    shown = ["fixed_rank1","energy_0.9_rank1","energy_0.5_rank1","correlated_rank2","radial_selected_rank2","aligned_rank1"]
    fig, axes = plt.subplots(2,3,figsize=(15,8),layout="constrained")
    for ax, case in zip(axes.flat, shown, strict=True):
        evaluation, grid, parameters = plot_data[case]
        values = np.sort(evaluation)
        indices = np.sort(len(values)-np.unique(np.geomspace(1,len(values),600).astype(int)))
        ax.semilogy(values[indices], (len(values)-indices)/len(values), color="0.7",label="Evaluation data",lw=3)
        for method,label,color,style in [
            ("current_chi_square","Current chi-square","#ba3d35","--"),
            ("scale_and_df_moment_fit","Fitted scale + df","#c88b13",":"),
        ]:
            a,nu=parameters[method]
            ax.semilogy(grid,chi2.sf(grid/a,nu),style,color=color,label=label,lw=2)
        ax.semilogy(grid,exact_sf(case,grid),color="#176b69",label="Exact example law",lw=1.6)
        ax.set(title=LABELS[case],xlabel="Test statistic",ylabel="Upper-tail probability",ylim=(1e-4,1.05),xlim=(0,grid[-1]))
        ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    fig.suptitle("Current, independently fitted, and exact distributions\nKnown Gaussian examples; no full-pipeline calibration claim",fontsize=14)
    fig.savefig(BASE/"distribution_comparison.png",dpi=170)
    fig.savefig(BASE/"distribution_comparison.pdf")
    plt.close(fig)
    table=pd.DataFrame(comparison_rows)
    print(table[table.alpha.eq(.01) & table.method.isin(["current_chi_square","scale_and_df_moment_fit","exact_conditional_law"])][["case","method","n_test","rejections","rate","true_rate_for_frozen_fit"]].to_string(index=False))


if __name__ == "__main__":
    main()
