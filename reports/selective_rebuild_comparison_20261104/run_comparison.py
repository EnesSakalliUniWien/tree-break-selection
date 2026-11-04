"""Bounded root-only replacement diagnostics; see protocol.md for the fixed scope."""

import hashlib
import itertools
import json
import platform
import sys
import time
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import chi, norm
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import (
    _geometry,
    make_diffusion_root_path,
)

HERE = Path(__file__).resolve().parent
SCENARIOS = {"null": (12, 0.0), "moderate": (6, 1.0), "strong": (6, 2.0)}


def root_groups(signature, n):
    groups = {i: (i,) for i in range(n)}
    for step, (a, b) in enumerate(signature):
        groups[n + step] = tuple(sorted(groups[a] + groups[b]))
    a, b = signature[-1]
    return tuple(sorted((groups[a], groups[b])))


def split_score(values, signature):
    groups = root_groups(signature, len(values))
    contrast = np.zeros(len(values))
    contrast[list(groups[0])] = 1 / len(groups[0])
    contrast[list(groups[1])] = -1 / len(groups[1])
    contrast /= np.linalg.norm(contrast)
    bss = float(np.sum((contrast @ values) ** 2))
    tss = float(np.sum((values - values.mean(axis=0)) ** 2))
    if tss <= 0:
        raise ValueError("undefined cluster index for constant observations")
    return bss / tss


def geometry_score(values):
    signature, _ = _geometry(values, tuple(f"S{i}" for i in range(len(values))))
    return split_score(values, signature)


def importance_summary(radii, selected, observed, dimension):
    # The negative normal tail has zero target mass and never reaches geometry.
    valid = selected & (radii > 0)
    if not valid.any():
        return dict(p_value=None, ess=0.0, mcse=None, accepted=0)
    r = radii[valid]
    log_target = chi.logpdf(r, dimension)
    log_proposal = np.logaddexp(log_target, norm.logpdf(r, loc=observed)) - np.log(2)
    logw = log_target - log_proposal
    weights = np.exp(logw - logsumexp(logw))
    upper = r >= observed
    estimate = float(weights @ upper)
    return dict(
        p_value=estimate,
        ess=float(1 / (weights @ weights)),
        mcse=float(np.sqrt(np.sum(weights**2 * (upper - estimate) ** 2))),
        accepted=int(valid.sum()),
    )


def chain_rank(observed, dimension, event, *, steps, seed):
    rng = np.random.default_rng(seed)
    anchor = int(rng.integers(steps + 1))
    states = np.full(steps + 1, np.nan)
    states[anchor] = observed
    moves = 0
    evaluations = 0
    for indices, source_delta in [
        (range(anchor - 1, -1, -1), 1),
        (range(anchor + 1, steps + 1), -1),
    ]:
        for index in indices:
            current = states[index + source_delta]
            proposed = current + rng.normal()
            states[index] = current
            if proposed <= 0:
                continue
            log_ratio = (dimension - 1) * np.log(proposed / current) - (
                proposed**2 - current**2
            ) / 2
            if np.log(rng.random()) > min(0.0, log_ratio):
                continue
            evaluations += 1
            if event(proposed):
                states[index] = proposed
                moves += 1
    assert np.isfinite(states).all() and states[anchor] == observed
    return dict(
        p_value=float(np.mean(states >= observed)),
        moves=moves,
        unique_states=int(len(np.unique(states))),
        geometry_evaluations=evaluations,
        anchor=anchor,
    ), states


def self_check():
    first = ((0, 1), (2, 3), (4, 5))
    second = ((2, 3), (0, 1), (4, 5))
    assert first != second and root_groups(first, 4) == root_groups(second, 4)
    x = np.array([[0, 1], [1, 2], [4, 2], [5, 3]], dtype=float)
    within = sum(np.sum((x[list(g)] - x[list(g)].mean(axis=0)) ** 2) for g in root_groups(first, 4))
    np.testing.assert_allclose(
        split_score(x, first), 1 - within / np.sum((x - x.mean(axis=0)) ** 2)
    )
    rng = np.random.default_rng(2026110411)
    r = np.where(
        rng.random(200000) < 0.5, np.sqrt(rng.chisquare(3, 200000)), rng.normal(6, 1, 200000)
    )
    result = importance_summary(r, r >= 5, 6, 3)
    exact = float(chi.sf(6, 3) / chi.sf(5, 3))
    assert abs(result["p_value"] - exact) < 5 * result["mcse"] + 1e-4
    # Enumerate the random-position rank law for a nontrivial reversible chain.
    pi = np.array([0.3, 0.7])
    transition = np.array([[0.8, 0.2], [0.06 / 0.7, 1 - 0.06 / 0.7]])
    masses = []
    for path in itertools.product([0, 1], repeat=4):
        probability = pi[path[0]] * np.prod([transition[a, b] for a, b in zip(path[:-1], path[1:])])
        for anchor in range(4):
            p = np.mean(np.asarray(path) >= path[anchor])
            masses.append((p, probability / 4))
    assert abs(sum(m for _, m in masses) - 1) < 1e-12
    for alpha in [0.25, 0.5, 0.75, 1.0]:
        assert sum(m for p, m in masses if p <= alpha) <= alpha + 1e-12
    stalled, states = chain_rank(6, 3, lambda r: False, steps=99, seed=11)
    assert stalled["p_value"] == 1 and stalled["moves"] == 0 and np.all(states == 6)
    return dict(
        importance_exact=exact,
        importance_estimate=result,
        finite_state_chain="passed",
        root_and_score="passed",
    )


def evaluate_case(scenario, replicate, shift, reference_scores):
    started = time.perf_counter()
    seed = 2026110401 + replicate
    x = np.random.default_rng(seed).normal(size=(9, 3))
    means = np.zeros_like(x)
    means[:4, 0], means[4:, 0] = -shift, shift
    x += means
    path = make_diffusion_root_path(x, sigma=1)
    observed_root = root_groups(path.observed_signature, len(x))
    arrays = {"input": x, "true_means": means}
    records = []
    base = dict(
        scenario=scenario,
        replicate=replicate,
        data_seed=seed,
        observed_radius=path.observed_radius,
        true_root_contrast=float(np.linalg.norm(path.e @ means)),
        root_sizes=[len(g) for g in observed_root],
    )

    def membership(radius):
        signature, _ = _geometry(path.offset + radius * path.slope, path.labels)
        full = signature == path.observed_signature
        root = root_groups(signature, len(x)) == observed_root
        assert not full or root
        return full, root

    assert membership(path.observed_radius) == (True, True)
    rng = np.random.default_rng(2026120401 + replicate)
    draws = np.sqrt(rng.chisquare(3, 499))
    accepted = np.array([membership(float(r)) for r in draws], dtype=bool)
    arrays.update(naive_radii=draws, naive_membership=accepted)
    for index, event in enumerate(["full_history", "root_split"]):
        kept = draws[accepted[:, index]]
        records.append(
            dict(
                base,
                method="conditional_rejection",
                event=event,
                p_value=float(
                    (1 + np.count_nonzero(kept >= path.observed_radius)) / (1 + len(kept))
                ),
                accepted=len(kept),
                minimum_p_value=1 / (1 + len(kept)),
                status="numerical_rank",
            )
        )
    rng = np.random.default_rng(2026130401 + replicate)
    draws = np.where(
        rng.random(999) < 0.5,
        np.sqrt(rng.chisquare(3, 999)),
        rng.normal(path.observed_radius, 1, 999),
    )
    accepted = np.array(
        [membership(float(r)) if r > 0 else (False, False) for r in draws], dtype=bool
    )
    arrays.update(importance_radii=draws, importance_membership=accepted)
    for index, event in enumerate(["full_history", "root_split"]):
        result = importance_summary(draws, accepted[:, index], path.observed_radius, 3)
        records.append(
            dict(base, method="conditional_importance", event=event, status="approximate", **result)
        )
    for index, event in enumerate(["full_history", "root_split"]):
        result, states = chain_rank(
            path.observed_radius,
            3,
            lambda r, i=index: membership(r)[i],
            steps=999,
            seed=2026140401 + replicate,
        )
        arrays[f"chain_{event}"] = states
        records.append(
            dict(base, method="conditional_chain", event=event, status="numerical_rank", **result)
        )
    score = split_score(x, path.observed_signature)
    records.append(
        dict(
            base,
            method="rebuild_known_covariance",
            event="one_gaussian_population",
            p_value=float(
                (1 + np.count_nonzero(reference_scores >= score)) / (1 + len(reference_scores))
            ),
            status="numerical_rank",
            score=score,
        )
    )
    covariance = np.cov(x, rowvar=False)
    factor = np.linalg.cholesky(covariance)
    rng = np.random.default_rng(2026150401 + replicate)
    fitted_scores = np.array(
        [geometry_score(rng.normal(size=x.shape) @ factor.T) for _ in range(499)]
    )
    arrays.update(fitted_covariance=covariance, fitted_null_scores=fitted_scores)
    records.append(
        dict(
            base,
            method="rebuild_fitted_covariance",
            event="one_gaussian_population",
            p_value=float(
                (1 + np.count_nonzero(fitted_scores >= score)) / (1 + len(fitted_scores))
            ),
            status="plugin_approximate",
            score=score,
        )
    )
    for record in records:
        record["case_seconds"] = time.perf_counter() - started
    return records, arrays


def main():
    checks = self_check()
    if "--self-check" in sys.argv:
        print(json.dumps(checks, indent=2))
        return
    output = HERE / "generated"
    output.mkdir(parents=True, exist_ok=True)
    sources = [
        *ROOT.glob("tree_break_selection/**/*.py"),
        Path(__file__),
        ROOT / "benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_replay.py",
        HERE / "protocol.md",
    ]
    manifest = dict(
        python=sys.version,
        platform=platform.platform(),
        versions={
            p: version(p) for p in ["numpy", "scipy", "pandas", "pydiffmap", "threadpoolctl"]
        },
        scenarios=SCENARIOS,
        source_sha256={
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources
        },
        self_check=checks,
    )
    (HERE / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    records, errors = [], []
    with threadpool_limits(limits=1), (output / "records.jsonl").open("x") as stream:
        rng = np.random.default_rng(2026160401)
        reference = np.array([geometry_score(rng.normal(size=(9, 3))) for _ in range(999)])
        np.save(output / "known_null_scores.npy", reference)
        print("known reference bank complete", flush=True)
        for scenario, (replicates, shift) in SCENARIOS.items():
            for replicate in range(replicates):
                try:
                    case_records, arrays = evaluate_case(scenario, replicate, shift, reference)
                except (ValueError, RuntimeError, np.linalg.LinAlgError) as error:
                    errors.append(dict(scenario=scenario, replicate=replicate, error=str(error)))
                    print("UNRESOLVED", errors[-1], flush=True)
                    continue
                np.savez_compressed(output / f"{scenario}_{replicate}.npz", **arrays)
                for record in case_records:
                    stream.write(json.dumps(record) + "\n")
                stream.flush()
                records.extend(case_records)
                print(
                    scenario,
                    replicate,
                    "complete",
                    round(case_records[0]["case_seconds"], 2),
                    flush=True,
                )
    table = pd.DataFrame(records)
    table.to_csv(HERE / "results.csv", index=False)
    summary = (
        table.groupby(["scenario", "method", "event"], sort=False)
        .agg(
            cases=("p_value", "size"),
            available=("p_value", "count"),
            reject_005=("p_value", lambda s: int((s <= 0.05).sum())),
            reject_001=("p_value", lambda s: int((s <= 0.01).sum())),
            median_p=("p_value", "median"),
        )
        .reset_index()
    )
    summary.to_csv(HERE / "summary.csv", index=False)
    (HERE / "errors.json").write_text(json.dumps(errors, indent=2) + "\n")
    for name, digest in manifest["source_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
