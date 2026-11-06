"""Matched native topology observations around the identity-self-distance fix.

Capture with ``--stage baseline`` only while the genuine old source is live;
the command refuses to overwrite a baseline. Run ``--stage compare`` after
the fix using the archived original input arrays, without legacy monkeypatches.
Root chi-square tails are descriptive unselected known-I diagnostics only.
"""

import argparse
import hashlib
import importlib.metadata
import itertools
import json
import platform
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
from scipy.stats import chi2

from benchmarks.diagnostics.calibration.diffusion_factorial import datasets, digest
from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import (
    _geometry,
    _topology,
)
from tree_break_selection.space_separation import diffusion
from tree_break_selection.tree.construction.build import _canonical_sample_order


def _operator(scaled_distances, epsilon):
    kernel = scaled_distances.copy()
    kernel.data = np.exp(-(kernel.data**2) / (4 * epsilon))
    weights = kernel.toarray()
    np.fill_diagonal(weights, 0)
    weights = np.maximum(weights, weights.T)
    degrees = weights.sum(axis=1)
    degrees[degrees == 0] = 1e-10
    return weights / np.sqrt(degrees[:, None] * degrees[None, :])


def _observe_geometry(values, labels):
    original_eigh = diffusion.eigh
    observed = {}

    def capture(matrix, *args, **kwargs):
        observed["operator"] = np.array(matrix, copy=True)
        return original_eigh(matrix, *args, **kwargs)

    # Observe the actual matrix sent to the native dense solver, unchanged.
    with patch.object(diffusion, "eigh", capture):
        signature, state = _geometry(values, labels)
    state["native_transition_operator"] = observed["operator"]
    return signature, state


def _record(values, labels, signature, state, case, replicate, seed, arm):
    n, p = values.shape
    groups = {i: [i] for i in range(n)}
    for step, (a, b) in enumerate(signature):
        groups[n + step] = groups[a] + groups[b]
    a, b = signature[-1]
    left, right = groups[a], groups[b]
    root_groups = sorted([sorted(str(labels[i]) for i in group) for group in (left, right)])
    statistic = float(
        np.sum((values[left].mean(0) - values[right].mean(0)) ** 2)
        / (1 / len(left) + 1 / len(right))
    )
    epsilon = state["outer_epsilon"]
    native = state["native_transition_operator"]
    reference = _operator(
        state["scaled_distances"], state["off_diagonal_positive_scaled_sq_median"]
    )
    return dict(
        case=case,
        replicate=replicate,
        seed=seed,
        arm=arm,
        n=n,
        p=p,
        data_sha256=hashlib.sha256(values.tobytes()).hexdigest(),
        signature=signature,
        clades=_topology(signature, n),
        root_groups=root_groups,
        n_left=len(left),
        n_right=len(right),
        known_identity_statistic=statistic,
        unselected_chi_p=float(chi2.sf(statistic, p)),
        epsilon=epsilon,
        positive_self_count=state["positive_self_scaled_distance_count"],
        zero_self_count=state["zero_self_scaled_distance_count"],
        actual_positive_median=state["actual_positive_scaled_sq_median"],
        off_diagonal_positive_median=state["off_diagonal_positive_scaled_sq_median"],
        native_vs_zero_self_reference_operator_max_abs=float(np.max(np.abs(native - reference))),
    )


def observe_dataset(values, labels, case, replicate, seed):
    """Observe original input and its floating radial reconstruction serially."""
    order = _canonical_sample_order(list(labels))
    labels = tuple(labels[i] for i in order)
    values = np.asarray(values, dtype=float)[order].copy()
    signature, state = _observe_geometry(values, labels)
    rows = [_record(values, labels, signature, state, case, replicate, seed, "original")]
    groups = {i: [i] for i in range(len(values))}
    for step, (a, b) in enumerate(signature):
        groups[len(values) + step] = groups[a] + groups[b]
    a, b = signature[-1]
    eta = np.zeros(len(values))
    eta[groups[a]] = 1 / len(groups[a])
    eta[groups[b]] = -1 / len(groups[b])
    e = eta / np.linalg.norm(eta)
    contrast = e @ values
    radius = np.linalg.norm(contrast)
    if not radius > 0:
        raise ValueError("Observed radial direction is undefined.")
    reconstructed = values - np.outer(e, contrast) + radius * np.outer(e, contrast / radius)
    rebuilt, rebuilt_state = _observe_geometry(reconstructed, labels)
    rows.append(
        _record(
            reconstructed, labels, rebuilt, rebuilt_state, case, replicate, seed, "reconstructed"
        )
    )
    return rows


def compare_records(before, after):
    """Keep merge-sequence ordering separate from clades and root membership."""
    keys = ("case", "replicate", "arm")
    old = {tuple(row[key] for key in keys): row for row in before}
    new = {tuple(row[key] for key in keys): row for row in after}
    if old.keys() != new.keys() or len(old) != len(before) or len(new) != len(after):
        raise ValueError("Comparison records must match uniquely by dataset and arm.")
    rows = []
    for key, previous in old.items():
        current = new[key]
        old_signature = tuple(map(tuple, previous["signature"]))
        new_signature = tuple(map(tuple, current["signature"]))
        old_clades = {tuple(clade) for clade in previous["clades"]}
        new_clades = {tuple(clade) for clade in current["clades"]}
        rows.append(
            dict(zip(keys, key))
            | dict(
                full_history_changed=old_signature != new_signature,
                clades_changed=old_clades != new_clades,
                root_partition_changed=previous["root_groups"] != current["root_groups"],
                clade_symmetric_difference=len(old_clades ^ new_clades),
                epsilon_before=previous["epsilon"],
                epsilon_after=current["epsilon"],
                positive_self_before=previous["positive_self_count"],
                positive_self_after=current["positive_self_count"],
                statistic_before=previous["known_identity_statistic"],
                statistic_after=current["known_identity_statistic"],
                unselected_chi_p_before=previous["unselected_chi_p"],
                unselected_chi_p_after=current["unselected_chi_p"],
            )
        )
    return pd.DataFrame(rows)


def _inputs(replicates):
    for case, replicate, seed, _, data, *_ in itertools.islice(datasets(replicates, 0), replicates):
        yield case, replicate, seed, data
    for replicate in range(3):
        seed = 20261105 + replicate
        yield (
            "small_gaussian",
            replicate,
            seed,
            pd.DataFrame(
                np.random.default_rng(seed).normal(size=(9, 1)),
                index=[f"S{i}" for i in range(9)],
                columns=["F0"],
            ),
        )
    for case, replicate, seed, _, data, *_ in itertools.islice(datasets(0, 1), 2):
        yield case, replicate, seed, data


def _provenance():
    sources = [
        Path(__file__),
        Path(diffusion.__file__),
        Path(datasets.__code__.co_filename),
        Path(_geometry.__code__.co_filename),
        Path("uv.lock"),
    ]
    return dict(
        python=platform.python_version(),
        versions={
            name: importlib.metadata.version(name)
            for name in ("numpy", "scipy", "pydiffmap", "scikit-learn", "pandas")
        },
        source_sha256={
            str(source.resolve()): hashlib.sha256(source.read_bytes()).hexdigest()
            for source in sources
        },
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", choices=("baseline", "compare"), default="compare")
    parser.add_argument("--replicates", type=int, default=50)
    parser.add_argument("--baseline-dir", type=Path)
    parser.add_argument(
        "--output", type=Path, default=Path("raw/assets/diffusion-self-distance-20261102")
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if args.stage == "baseline":
        if (args.output / "baseline_runs.json").exists():
            parser.error("Refusing to overwrite a genuine pre-fix baseline.")
        arrays, metadata, rows = {}, [], []
        for k, (case, replicate, seed, data) in enumerate(_inputs(args.replicates)):
            order = _canonical_sample_order(data.index.tolist())
            values = data.to_numpy(dtype=float)[order].copy()
            labels = tuple(data.index[i] for i in order)
            key = f"case_{k}"
            arrays[key] = values
            metadata.append(
                dict(
                    key=key,
                    case=case,
                    replicate=replicate,
                    seed=seed,
                    labels=labels,
                    columns=list(map(str, data.columns)),
                    data_sha256=hashlib.sha256(values.tobytes()).hexdigest(),
                )
            )
            rows.extend(observe_dataset(values, labels, case, replicate, seed))
        np.savez_compressed(args.output / "baseline_inputs.npz", **arrays)
        (args.output / "baseline_inputs.json").write_text(json.dumps(metadata, indent=2) + "\n")
        (args.output / "baseline_runs.json").write_text(json.dumps(rows, indent=2) + "\n")
        (args.output / "pre_fix_diffusion.py").write_bytes(Path(diffusion.__file__).read_bytes())
        (args.output / "baseline_manifest.json").write_text(
            json.dumps(_provenance(), indent=2) + "\n"
        )
        return
    if (args.output / "post_fix_runs.json").exists():
        parser.error(
            "Refusing to overwrite raw comparison evidence; choose a new --output directory."
        )
    baseline_dir = args.baseline_dir or args.output
    baseline = json.loads((baseline_dir / "baseline_runs.json").read_text())
    metadata = json.loads((baseline_dir / "baseline_inputs.json").read_text())
    rows = []
    with np.load(baseline_dir / "baseline_inputs.npz") as arrays:
        for item in metadata:
            values = arrays[item["key"]]
            assert hashlib.sha256(values.tobytes()).hexdigest() == item["data_sha256"]
            rows.extend(
                observe_dataset(
                    values, item["labels"], item["case"], item["replicate"], item["seed"]
                )
            )
            print(f"compared {item['case']} replicate={item['replicate']}", flush=True)
    comparison = compare_records(baseline, rows)
    (args.output / "post_fix_runs.json").write_text(json.dumps(rows, indent=2) + "\n")
    comparison.to_csv(args.output / "matched_topology_changes.csv", index=False)
    stability = []
    for stage, records in (("before", baseline), ("after", rows)):
        original = [dict(row, arm="paired") for row in records if row["arm"] == "original"]
        reconstructed = [
            dict(row, arm="paired") for row in records if row["arm"] == "reconstructed"
        ]
        frame = compare_records(original, reconstructed)
        frame["stage"] = stage
        stability.append(frame)
    pd.concat(stability, ignore_index=True).to_csv(
        args.output / "reconstruction_stability.csv", index=False
    )
    summary = (
        comparison.groupby(["case", "arm"])
        .agg(
            datasets=("replicate", "size"),
            history_changes=("full_history_changed", "sum"),
            clade_changes=("clades_changed", "sum"),
            root_changes=("root_partition_changed", "sum"),
            positive_self_before=("positive_self_before", "sum"),
            positive_self_after=("positive_self_after", "sum"),
        )
        .reset_index()
    )
    rejections = (
        comparison.assign(
            known_i_unselected_chi_reject_before=comparison.unselected_chi_p_before < 0.01,
            known_i_unselected_chi_reject_after=comparison.unselected_chi_p_after < 0.01,
        )
        .groupby(["case", "arm"])[
            ["known_i_unselected_chi_reject_before", "known_i_unselected_chi_reject_after"]
        ]
        .sum()
        .reset_index()
    )
    summary = summary.merge(rejections, on=["case", "arm"], validate="one_to_one")
    summary.to_csv(args.output / "comparison_summary.csv", index=False)
    manifest = _provenance() | dict(
        stage="post_fix_comparison",
        datasets=len(metadata),
        rows=len(rows),
        baseline_source_sha256=hashlib.sha256(
            (baseline_dir / "pre_fix_diffusion.py").read_bytes()
        ).hexdigest(),
        input_arrays_sha256=hashlib.sha256(
            (baseline_dir / "baseline_inputs.npz").read_bytes()
        ).hexdigest(),
        baseline_directory=str(baseline_dir),
        model="full-feature known-I root statistic; no NNLS/PCA/fitted inflation; unselected chi tails are descriptive only",
        operator_reference="floating zero-self/positive-offdiagonal epsilon reduction with native row-dependent bandwidths; not a rigorous exact-real certificate",
        baseline_operator_metric="not captured in original pre-fix rows; no legacy monkeypatch used to fill it",
        source_changes="baseline is archived genuine live pre-fix observations; inputs replayed under live fixed source",
    )
    prior_manifest = Path("raw/assets/diffusion-factorial-20261102/manifest.json")
    manifest["prior_factorial_manifest_sha256"] = hashlib.sha256(
        prior_manifest.read_bytes()
    ).hexdigest()
    prior = pd.read_csv(prior_manifest.with_name("runs.csv"))
    prior = (
        prior[(prior["case"] == "iid_continuous") & (prior.diffusion == 1)]
        .groupby("replicate")
        .first()
    )
    matches = 0
    for row in baseline:
        if row["case"] == "iid_continuous" and row["arm"] == "original":
            assert (
                digest([[float(a), float(b)] for a, b in row["signature"]])
                == prior.loc[row["replicate"], "topology_hash"]
            )
            matches += 1
    manifest["genuine_baseline_histories_matching_prior_factorial"] = matches
    manifest["raw_preservation_note"] = (
        "First postfix files were overwritten once while adding direct operator observation, before overwrite protection was implemented. Genuine pre-fix baseline is unchanged. Existing root-level postfix outputs are retained; this version is a separate additional run."
    )
    (args.output / "comparison_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
