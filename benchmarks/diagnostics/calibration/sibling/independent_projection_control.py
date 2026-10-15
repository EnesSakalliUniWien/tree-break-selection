"""Matched, independent-data controls for PCA basis and MP rank selection."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import shlex
import subprocess
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import chi2

from benchmarks.diagnostics.calibration.sibling.known_branch_skew_diagnostic import (
    INTERNAL_CHILDREN,
    LEAF_NAMES,
    PATTERNS,
    PROFILES,
    KnownBranchConfig,
    covariance_matrix,
    expected_node_contrasts,
    population_means,
)
from benchmarks.diagnostics.calibration.statistics.statistic_distribution_shape_panel import (
    _expected_chi_square_mixture_skew,
)
from benchmarks.validation.contracts.report_contract import wilson_interval
from tree_break_selection.hierarchy_analysis.decomposition.backends.eigen.decomposition import (
    eigendecompose_covariance,
)
from tree_break_selection.hierarchy_analysis.decomposition.backends.eigen.projection import (
    build_pca_projection,
)
from tree_break_selection.hierarchy_analysis.statistics.contrast_covariance import (
    build_null_whitened_tangent_matrix,
    compute_whitened_wald_contrast,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald.projected_wald_kernel import (
    run_projected_wald_kernel,
)
from tree_break_selection.hierarchy_analysis.statistics.projection.projection_dimension_estimation.projection_dimension_estimators import (
    estimate_marchenko_pastur_dimension,
)
from tree_break_selection.tree.feature_space import continuous_feature_space_from_columns

ROOT = Path(__file__).resolve().parents[4]
BASIS_SOURCES = ("evaluation", "independent", "fixed_coordinates")
GROUP = ["regime", "pattern", "covariance_profile", "focal_node"]
TRAINING_NAMESPACE = 17021983


@dataclass(frozen=True)
class ControlConfig:
    replicates: int = 1000
    base_seed: int = 20260928
    n_per_group: int = 32
    dimension: int = 8
    alpha: float = 0.01
    patterns: tuple[str, ...] = PATTERNS
    covariance_profiles: tuple[str, ...] = PROFILES
    stress_threshold_multiplier: int = 64

    def __post_init__(self):
        KnownBranchConfig(
            output_dir=Path("."),
            **{k: v for k, v in asdict(self).items() if k != "stress_threshold_multiplier"},
        )
        if self.stress_threshold_multiplier <= 1:
            raise ValueError("Stress multiplier must exceed one.")


def data_hash(data):
    return hashlib.sha256(np.ascontiguousarray(data).tobytes()).hexdigest()


def draw_data(seed, means, sigma, n_per_group):
    rng = np.random.default_rng(seed)
    return np.vstack(
        [rng.multivariate_normal(means[leaf], sigma, size=n_per_group) for leaf in LEAF_NAMES]
    )


def node_data(data):
    a0, a1, b0, b1 = np.split(data, 4)
    return {
        "A0": a0,
        "A1": a1,
        "B0": b0,
        "B1": b1,
        "A": np.vstack((a0, a1)),
        "B": np.vstack((b0, b1)),
        "R": data,
    }


def spectrum(data, sigma, feature_space):
    tangent = build_null_whitened_tangent_matrix(
        data,
        data.mean(axis=0),
        feature_space=feature_space,
        continuous_covariance_by_block={"continuous": sigma},
        ridge=0,
    )
    eigen = eigendecompose_covariance(tangent, compute_eigenvectors=True)
    if eigen is None:
        raise ValueError("Zero-rank tangent matrix is unsupported.")
    return eigen


def mp_rank(eigen, n_rows, multiplier=1):
    return estimate_marchenko_pastur_dimension(
        eigen.eigenvalues,
        n_samples=n_rows,
        n_features=eigen.active_feature_count,
        effective_independent_rows=n_rows,
        mp_threshold_rows=n_rows * multiplier,
        minimum_projection_dimension=1,
    )


def evaluate_projection(z, eigen, rank):
    """Fixed coordinate prefixes bypass PCA ties; PCA retains complete groups."""
    if eigen is None:
        statistic = float(z[:rank] @ z[:rank])
        return dict(statistic=statistic, p_value=float(chi2.sf(statistic, rank)), actual_rank=rank)
    projection, eigenvalues = build_pca_projection(
        eigen, projection_dimension=rank, n_features_total=len(z)
    )
    result = run_projected_wald_kernel(
        z, spectral_k=projection.shape[0], pca_projection=projection, pca_eigenvalues=eigenvalues
    )
    return dict(
        statistic=result.statistic, p_value=result.p_value, actual_rank=result.projection_dimension
    )


def run_control(config):
    rows, datasets = [], []
    seeds = iter(
        np.random.SeedSequence(config.base_seed).spawn(
            len(config.patterns) * len(config.covariance_profiles) * config.replicates
        )
    )
    feature_space = continuous_feature_space_from_columns(
        tuple(f"x{i}" for i in range(config.dimension))
    )
    for profile in config.covariance_profiles:
        sigma = covariance_matrix(profile, config.dimension)
        for pattern in config.patterns:
            means = population_means(pattern, config.dimension)
            truths = expected_node_contrasts(means)
            for replicate in range(config.replicates):
                evaluation_seed = next(seeds)
                # Original stream is unchanged. Training keys cannot collide with its one-level keys.
                basis_seed = np.random.SeedSequence(
                    config.base_seed, spawn_key=(TRAINING_NAMESPACE, 1, *evaluation_seed.spawn_key)
                )
                rank_seed = np.random.SeedSequence(
                    config.base_seed, spawn_key=(TRAINING_NAMESPACE, 2, *evaluation_seed.spawn_key)
                )
                sources = {
                    name: draw_data(seed, means, sigma, config.n_per_group)
                    for name, seed in (
                        ("evaluation", evaluation_seed),
                        ("basis_training", basis_seed),
                        ("rank_training", rank_seed),
                    )
                }
                common = dict(
                    pattern=pattern,
                    covariance_profile=profile,
                    replicate=replicate,
                    dataset_sha256=data_hash(sources["evaluation"]),
                )
                datasets.append(
                    {
                        **common,
                        **{f"{name}_sha256": data_hash(data) for name, data in sources.items()},
                        "evaluation_spawn_key": json.dumps(evaluation_seed.spawn_key),
                        "basis_training_spawn_key": json.dumps(basis_seed.spawn_key),
                        "rank_training_spawn_key": json.dumps(rank_seed.spawn_key),
                    }
                )
                nodes = {name: node_data(data) for name, data in sources.items()}
                for node, (left, right) in INTERNAL_CHILDREN.items():
                    contrast, is_null = truths[node]
                    base = {
                        **common,
                        "focal_node": node,
                        "population_null": is_null,
                        "population_contrast": json.dumps(contrast.tolist()),
                        "n_rows": len(nodes["evaluation"][node]),
                    }
                    failure = ""
                    try:
                        nl, nr = nodes["evaluation"][left], nodes["evaluation"][right]
                        z = compute_whitened_wald_contrast(
                            nl.mean(axis=0),
                            nr.mean(axis=0),
                            len(nl),
                            len(nr),
                            comparison="sibling",
                            feature_space=feature_space,
                            continuous_covariance_by_block={"continuous": sigma},
                            ridge=0,
                        )
                        eigens = {
                            name: spectrum(data[node], sigma, feature_space)
                            for name, data in nodes.items()
                        }
                    except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
                        failure = str(exc)
                    regimes = [("canonical_mp", 1)]
                    if pattern == "global_null":
                        regimes.append(
                            ("threshold_rows_stress", config.stress_threshold_multiplier)
                        )
                    for regime, multiplier in regimes:
                        rank_sources = ["evaluation_mp", "independent_mp"]
                        if regime == "canonical_mp":
                            rank_sources += [f"fixed_{k}" for k in sorted({1, 2, config.dimension})]
                        for basis_source in BASIS_SOURCES:
                            for rank_source in rank_sources:
                                row = {
                                    **base,
                                    "regime": regime,
                                    "basis_source": basis_source,
                                    "rank_source": rank_source,
                                    "mp_threshold_rows": base["n_rows"] * multiplier,
                                    "status": "runtime_error",
                                    "failure_reason": failure,
                                    "requested_rank": math.nan,
                                    "actual_rank": math.nan,
                                    "raw_mp_signal_count": math.nan,
                                    "statistic": math.nan,
                                    "p_value": math.nan,
                                    "conditional_chi_square_guarantee": False,
                                }
                                if not failure:
                                    try:
                                        if rank_source.startswith("fixed_"):
                                            rank = int(rank_source.removeprefix("fixed_"))
                                            raw_count = math.nan
                                        else:
                                            rank_eigen = eigens[
                                                "evaluation"
                                                if rank_source == "evaluation_mp"
                                                else "rank_training"
                                            ]
                                            estimate = mp_rank(
                                                rank_eigen, base["n_rows"], multiplier
                                            )
                                            rank, raw_count = (
                                                estimate.test_projection_dimension,
                                                estimate.raw_mp_signal_count,
                                            )
                                        eigen = (
                                            None
                                            if basis_source == "fixed_coordinates"
                                            else eigens[
                                                "evaluation"
                                                if basis_source == "evaluation"
                                                else "basis_training"
                                            ]
                                        )
                                        result = evaluate_projection(z, eigen, rank)
                                        row.update(
                                            result,
                                            requested_rank=rank,
                                            raw_mp_signal_count=raw_count,
                                            status="ok",
                                        )
                                        # Full predetermined rank is invariant to the selected basis.
                                        guarantee = (
                                            basis_source != "evaluation"
                                            and rank_source != "evaluation_mp"
                                        ) or (
                                            rank_source == f"fixed_{config.dimension}"
                                            and result["actual_rank"] == config.dimension
                                        )
                                        row["conditional_chi_square_guarantee"] = bool(
                                            is_null and guarantee
                                        )
                                    except (
                                        ValueError,
                                        np.linalg.LinAlgError,
                                        FloatingPointError,
                                    ) as exc:
                                        row["failure_reason"] = str(exc)
                                rows.append(row)
            print(f"completed {profile}/{pattern}: {config.replicates} datasets", flush=True)
    return pd.DataFrame(rows), pd.DataFrame(datasets)


def summarize(rows, alpha):
    results = []
    for key, frame in rows.groupby(GROUP + ["basis_source", "rank_source"], sort=True):
        ok = frame[frame.status == "ok"]
        rejects = int((ok.p_value <= alpha).sum())
        lo, hi = wilson_interval(rejects, len(ok)) if len(ok) else (math.nan, math.nan)
        ranks = ok.actual_rank.value_counts().sort_index()
        results.append(
            {
                **dict(zip(GROUP + ["basis_source", "rank_source"], key)),
                "population_null": bool(frame.population_null.iloc[0]),
                "n_total": len(frame),
                "n_ok": len(ok),
                "n_unavailable": len(frame) - len(ok),
                "status_counts": json.dumps(frame.status.value_counts().to_dict()),
                "rejections": rejects,
                "rejection_rate": rejects / len(ok) if len(ok) else math.nan,
                "wilson95_low": lo,
                "wilson95_high": hi,
                "rank_counts": json.dumps({str(int(k)): int(v) for k, v in ranks.items()}),
                "raw_mp_count_counts": json.dumps(
                    {
                        str(int(k)): int(v)
                        for k, v in ok.raw_mp_signal_count.dropna()
                        .value_counts()
                        .sort_index()
                        .items()
                    }
                ),
                "rank_varies": len(ranks) > 1,
                "mean_statistic_per_rank": (ok.statistic / ok.actual_rank).mean(),
                "null_statistic_skew": ok.statistic.skew()
                if bool(frame.population_null.iloc[0])
                else math.nan,
                "chi_square_df_mixture_reference_skew": _expected_chi_square_mixture_skew(
                    ok.actual_rank
                )
                if bool(frame.population_null.iloc[0]) and len(ok)
                else math.nan,
                "p_value_q01": ok.p_value.quantile(0.01),
                "p_value_q05": ok.p_value.quantile(0.05),
                "p_value_median": ok.p_value.median(),
                "conditional_chi_square_guarantee": bool(
                    frame.conditional_chi_square_guarantee.all()
                ),
            }
        )
    return pd.DataFrame(results)


def paired_outcomes(rows, alpha):
    summaries, outcomes = [], []
    keys = GROUP + ["replicate", "dataset_sha256"]
    comparisons = []
    for rank in sorted(rows.rank_source.unique()):
        comparisons += [
            ("basis_evaluation_minus_independent", "evaluation", rank, "independent", rank),
            ("basis_evaluation_minus_fixed", "evaluation", rank, "fixed_coordinates", rank),
        ]
    for basis in BASIS_SOURCES:
        comparisons.append(
            ("rank_evaluation_minus_independent", basis, "evaluation_mp", basis, "independent_mp")
        )
    for name, ba, ra, bb, rb in comparisons:
        a = rows[(rows.basis_source == ba) & (rows.rank_source == ra)]
        b = rows[(rows.basis_source == bb) & (rows.rank_source == rb)]
        joined = a.merge(b, on=keys, suffixes=("_a", "_b"), validate="one_to_one")
        if joined.empty:
            continue
        both = (joined.status_a == "ok") & (joined.status_b == "ok")
        joined["pair_status"] = np.where(both, "ok", "unavailable")
        joined["a_rejects"] = np.where(both, joined.p_value_a <= alpha, np.nan)
        joined["b_rejects"] = np.where(both, joined.p_value_b <= alpha, np.nan)
        joined["rejection_difference"] = joined.a_rejects - joined.b_rejects
        joined["statistic_difference"] = joined.statistic_a - joined.statistic_b
        joined["comparison"] = name
        cols = keys + [
            "comparison",
            "basis_source_a",
            "rank_source_a",
            "basis_source_b",
            "rank_source_b",
            "pair_status",
            "a_rejects",
            "b_rejects",
            "rejection_difference",
            "statistic_difference",
            "actual_rank_a",
            "actual_rank_b",
        ]
        outcomes.append(joined[cols])
        for key, group in joined.groupby(GROUP, sort=True):
            ok = group[group.pair_status == "ok"]
            diff = ok.rejection_difference
            se = diff.std(ddof=1) / np.sqrt(len(diff)) if len(diff) > 1 else math.nan
            summaries.append(
                {
                    **dict(zip(GROUP, key)),
                    "comparison": name,
                    "basis_source_a": ba,
                    "rank_source_a": ra,
                    "basis_source_b": bb,
                    "rank_source_b": rb,
                    "n_total": len(group),
                    "n_pairs": len(ok),
                    "n_unavailable": len(group) - len(ok),
                    "both_reject": int(((ok.a_rejects == 1) & (ok.b_rejects == 1)).sum()),
                    "a_only_reject": int(((ok.a_rejects == 1) & (ok.b_rejects == 0)).sum()),
                    "b_only_reject": int(((ok.a_rejects == 0) & (ok.b_rejects == 1)).sum()),
                    "neither_reject": int(((ok.a_rejects == 0) & (ok.b_rejects == 0)).sum()),
                    "rejection_difference": diff.mean(),
                    "paired_normal95_low": diff.mean() - 1.96 * se,
                    "paired_normal95_high": diff.mean() + 1.96 * se,
                    "actual_rank_disagreement": int((ok.actual_rank_a != ok.actual_rank_b).sum()),
                    "mean_statistic_difference": ok.statistic_difference.mean(),
                }
            )
    return pd.concat(outcomes, ignore_index=True), pd.DataFrame(summaries)


def source_hashes(protocol):
    paths = {Path(__file__), protocol}
    # Every imported project module that contributes to this run is captured,
    # including transitive helpers; dirty sources are hashed as actually used.
    for module in list(sys.modules.values()):
        path = getattr(module, "__file__", None)
        if path:
            p = Path(path).resolve()
            if p.is_relative_to(ROOT) and p.suffix == ".py" and ".venv" not in p.parts:
                paths.add(p)
    return {
        str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)
    }


def freeze(output, config, protocol, command):
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Refusing populated output: {output}")
    output.mkdir(parents=True, exist_ok=True)
    config_json = json.loads(json.dumps(asdict(config)))
    manifest = dict(
        schema="independent_projection_control/v1",
        frozen_utc=datetime.now(timezone.utc).isoformat(),
        config=config_json,
        protocol=str(protocol.relative_to(ROOT)),
        source_sha256=source_hashes(protocol),
        random_stream_policy="evaluation: SeedSequence(base_seed).spawn(profile_count*pattern_count*replicates); basis/rank training: spawn_key=(17021983,1/2,evaluation_key)",
        runtime={
            "python": sys.version,
            "platform": platform.platform(),
            **{p: version(p) for p in ("numpy", "scipy", "pandas")},
        },
        git_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        git_status=subprocess.check_output(["git", "status", "--short"], text=True).splitlines(),
        freeze_command=command,
    )
    (output / "frozen_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def verify_frozen(output, config, protocol):
    manifest = json.loads((output / "frozen_manifest.json").read_text())
    if manifest["config"] != json.loads(json.dumps(asdict(config))):
        raise ValueError("Configuration differs from frozen protocol.")
    for name, expected in manifest["source_sha256"].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Source changed after freeze: {name}")
    if str(protocol.relative_to(ROOT)) != manifest["protocol"]:
        raise ValueError("Protocol path differs from freeze.")
    if set(p.name for p in output.iterdir()) != {"frozen_manifest.json"}:
        raise FileExistsError("Run output already contains results; refusing overwrite.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("freeze", "run"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--base-seed", type=int, default=20260928)
    args = parser.parse_args()
    config = ControlConfig(replicates=args.replicates, base_seed=args.base_seed)
    output, protocol = args.output_dir.resolve(), args.protocol.resolve()
    command = shlex.join(
        [
            sys.executable,
            "-m",
            "benchmarks.diagnostics.calibration.sibling.independent_projection_control",
            *sys.argv[1:],
        ]
    )
    if args.action == "freeze":
        freeze(output, config, protocol, command)
        print(output / "frozen_manifest.json")
        return
    verify_frozen(output, config, protocol)
    started = datetime.now(timezone.utc).isoformat()
    rows, datasets = run_control(config)
    summary = summarize(rows, config.alpha)
    paired, paired_summary = paired_outcomes(rows, config.alpha)
    rows.to_csv(output / "control_rows.csv.gz", index=False)
    datasets.to_csv(output / "datasets.csv", index=False)
    summary.to_csv(output / "summary.csv", index=False)
    paired.to_csv(output / "paired_outcomes.csv.gz", index=False)
    paired_summary.to_csv(output / "paired_summary.csv", index=False)
    frozen = json.loads((output / "frozen_manifest.json").read_text())
    changed_sources = [
        name
        for name, expected in frozen["source_sha256"].items()
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != expected
    ]
    (output / "completion.json").write_text(
        json.dumps(
            dict(
                started_utc=started,
                completed_utc=datetime.now(timezone.utc).isoformat(),
                run_command=command,
                dataset_count=len(datasets),
                row_count=len(rows),
                source_hashes_unchanged=not changed_sources,
                changed_sources=changed_sources,
                status_counts=rows.status.value_counts().to_dict(),
                output_sha256={
                    p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in output.iterdir()
                    if p.is_file()
                },
            ),
            indent=2,
        )
        + "\n"
    )
    if changed_sources:
        raise RuntimeError(
            f"Sources changed during evaluation; outputs retained but unverified: {changed_sources}"
        )
    print(
        json.dumps(
            {
                "datasets": len(datasets),
                "rows": len(rows),
                "status": rows.status.value_counts().to_dict(),
            }
        )
    )


if __name__ == "__main__":
    main()
