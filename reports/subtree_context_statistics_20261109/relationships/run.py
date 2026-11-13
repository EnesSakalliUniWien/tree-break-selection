"""Join saved subspace coverage to raw covariance and native spectral magnitudes."""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from benchmarks.diagnostics.calibration.selected.hierarchy import (
    subtree_spectral_alignment as align,
)
from reports.subtree_spectral_alignment_20261109.alternatives.run import REGIMES, dataset
from tree_break_selection.eigenvalue_groups import tied_eigenvalue_neighbors

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PRIOR = ROOT / "reports/subtree_spectral_alignment_20261109"
ALTERNATIVES = PRIOR / "alternatives"
NATIVE = PRIOR / "parent_child"
TRACES = sorted(NATIVE.glob("trace_*.json.gz"))
INPUTS = (
    ALTERNATIVES / "rows.csv.gz",
    ALTERNATIVES / "summary.csv",
    ALTERNATIVES / "evidence.json",
    ALTERNATIVES / "run.py",
    ALTERNATIVES / "source_during_run.py.txt",
    ALTERNATIVES / "source_change_review.json",
    NATIVE / "rows.jsonl.gz",
    NATIVE / "manifest.json",
    NATIVE / "production_span_summary.json",
    *TRACES,
    Path(__file__),
    Path(align.__file__),
    Path(tied_eigenvalue_neighbors.__code__.co_filename),
)


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def inputs_hashes():
    return {str(path.relative_to(ROOT)): digest(path) for path in INPUTS}


def score_magnitudes():
    rows = pd.read_csv(ALTERNATIVES / "rows.csv.gz")
    assert len(rows) == 5120 and rows.status.eq("ok").all()
    enriched, paired = [], []
    spectra = {}
    for (regime, rep), original in rows.groupby(["regime", "replicate"], sort=False):
        x, truth, seed = dataset(regime, rep)
        data_hash = hashlib.sha256(x.tobytes()).hexdigest()
        assert original.seed.eq(seed).all() and original.data_sha256.eq(data_hash).all()
        spectrum = align.covariance_spectrum(x)
        values = spectrum.eigenvalues
        trace = float(values.sum())
        shares = values / trace
        spectra[(regime, rep)] = (data_hash, seed, values)
        group = original.copy()
        subspaces = group.blocks.map(ast.literal_eval)
        assert set(block for blocks in subspaces for block in blocks) == set(spectrum.blocks)
        group["covariance_trace"] = trace
        group["covariance_eigenvalue_sum"] = [
            float(sum(values[a:b].sum() for a, b in blocks)) if blocks else trace
            for blocks in subspaces
        ]
        group["covariance_trace_share"] = group.covariance_eigenvalue_sum / trace
        assert group.covariance_trace_share.between(0, 1 + 1e-12).all()
        enriched.append(group)
        if truth is None:
            continue
        full = group.loc[group.is_full_coordinates].iloc[0]
        leading = group.loc[group.is_leading].iloc[0]
        selected = group.loc[group.selected_common_target].iloc[0]
        candidates = group.loc[~group.is_full_coordinates]
        for factor, split, ami in (
            ("A", "split_a", "factor_a_ami"),
            ("B", "split_b", "contrast_b_ami"),
        ):
            covered = candidates.loc[candidates[split].eq(True)]
            ranked = covered.sort_values(["covariance_trace_share", "rank", "blocks"])
            minimum = ranked.iloc[0] if len(ranked) else None
            fit_covered = covered.loc[covered.shortlisted].sort_values("target_rmse")
            paired.append(
                {
                    "regime": regime,
                    "replicate": rep,
                    "seed": seed,
                    "data_hash": data_hash,
                    "factor": factor,
                    "full_split": full[split],
                    "leading_split": leading[split],
                    "selected_split": selected[split],
                    "full_ami": full[ami],
                    "leading_ami": leading[ami],
                    "selected_ami": selected[ami],
                    "family_max_ami_oracle": candidates[ami].max(),
                    "family_split": bool(len(ranked)),
                    "minimum_covering_share": minimum.covariance_trace_share
                    if minimum is not None
                    else np.nan,
                    "minimum_covering_blocks": minimum.blocks if minimum is not None else "",
                    "minimum_covering_rank": minimum["rank"] if minimum is not None else np.nan,
                    "minimum_covering_ami": minimum[ami] if minimum is not None else np.nan,
                    "leading_share": leading.covariance_trace_share,
                    "pc1_share": shares[0],
                    "pc2_share": shares[1],
                    "common_target_full_rmse": full.target_rmse,
                    "common_target_selected_rmse": selected.target_rmse,
                    "common_target_best_covering_rmse": fit_covered.target_rmse.iloc[0]
                    if len(fit_covered)
                    else np.nan,
                }
            )
    assert len(spectra) == 80 and len(paired) == 120
    return pd.concat(enriched, ignore_index=True), pd.DataFrame(paired), spectra


def native_magnitudes(spectra):
    with gzip.open(NATIVE / "rows.jsonl.gz", "rt") as stream:
        records = [json.loads(line) for line in stream]
    assert len(records) == 160 and len(TRACES) == 8
    rows, traced = [], []
    for item in records:
        case, rep, operator = (item[key] for key in ("case", "replicate", "operator"))
        data_hash, seed, values = spectra[(case, rep)]
        assert item["data_hash"] == data_hash and item["seed"] == seed
        assert item["status"] == "ok" and item["result"]["native"]["status"] == "ok"
        result = item["result"]
        raw = np.asarray(result["covariance"]["eigenvalues"])
        np.testing.assert_allclose(raw, values, rtol=1e-12, atol=1e-12)
        native = result["native"]
        root = result["native_nodes"][0]
        context = result["native_context"][root]
        whitened = np.asarray(context["principal_component_eigenvalues_by_node"])
        d = result["diffusion"]
        start, end = d["window"]
        spectrum = np.asarray(d["eigenvalues"])
        coefficients = np.maximum(spectrum[start:end], 0) ** 6
        total = float(coefficients.sum())
        first_block = next(
            (a, b) for item in d["whole_blocks"] for a, b in (item["block"],) if a <= start < b
        )
        first_block_count = min(end, first_block[1]) - start
        row = {
            "regime": case,
            "replicate": rep,
            "seed": seed,
            "data_hash": data_hash,
            "operator": operator,
            "root_p": native["root_sibling_p"],
            "root_reject": native["root_sibling_reject"],
            "native_ami": native["ami"],
            "native_clusters": native["found_clusters"],
            "root_first_child": repr(sorted(result["children"][0]["rows"])),
            "root_split": repr(sorted(sorted(child["rows"]) for child in result["children"])),
            "covariance_trace": float(raw.sum()),
            "raw_pc1_share": float(raw[0] / raw.sum()),
            "between_covariance_share": result["covariance"]["between_trace"] / raw.sum(),
            "diffusion_first_retained_lambda": spectrum[start],
            "diffusion_first_retained_coefficient_share": float(coefficients[0] / total),
            "diffusion_first_retained_block_rank": first_block_count,
            "diffusion_first_retained_block_coefficient_share": float(
                coefficients[:first_block_count].sum() / total
            ),
            "diffusion_coefficient_sum": total,
            "native_whitened_min_eig": float(whitened.min()),
            "native_whitened_max_eig": float(whitened.max()),
            "native_test_dimension": context["test_projection_dimensions_by_node"],
            "native_raw_mp_count": context["raw_mp_signal_counts_by_node"],
        }
        rows.append(row)
        if rep != 0:
            continue
        with gzip.open(NATIVE / f"trace_{case}_{operator}.json.gz", "rt") as stream:
            evidence = json.load(stream)
        spec = align.reversible_diffusion_spectrum(evidence["matrix_evidence"]["parent_weights"])
        np.testing.assert_allclose(spectrum, spec.eigenvalues, rtol=1e-12, atol=1e-12)
        basis = spec.symmetric_vectors if operator == "symmetric" else spec.right_vectors
        coordinates = basis[:, start:end] * np.sqrt(coefficients)
        component_energy = len(coordinates) * np.sum(
            (coordinates - coordinates.mean(axis=0)) ** 2, axis=0
        )
        squared_pair_energy = float(
            np.sum(np.asarray(evidence["baseline_diffusion_distance"]) ** 2)
        )
        np.testing.assert_allclose(
            component_energy.sum(), squared_pair_energy, rtol=1e-10, atol=1e-12
        )
        statistic = next(record for record in evidence["raw_records"] if record["parent"] == root)
        annotation = next(record for record in evidence["annotations"] if record["node_id"] == root)
        traced.append(
            row
            | {
                "diffusion_actual_first_mode_pair_energy_share": float(
                    component_energy[0] / component_energy.sum()
                ),
                "diffusion_actual_first_block_pair_energy_share": float(
                    component_energy[:first_block_count].sum() / component_energy.sum()
                ),
                "actual_pair_energy_total": squared_pair_energy,
                "sibling_raw_statistic": statistic["stat"],
                "sibling_reference_scale": statistic["reference_scale"],
                "sibling_degrees_of_freedom": statistic["degrees_of_freedom"],
                "sibling_applied_inflation": annotation["Sibling_Inflation_Factor"],
                "sibling_branch_variance_multiplier": 1
                + statistic["branch_length_sum"] / (2 * native["mean_positive_branch"]),
            }
        )
    native = pd.DataFrame(rows)
    left = native.loc[native.operator == "symmetric"].drop(columns="operator")
    right = native.loc[native.operator == "markov"].drop(columns="operator")
    paired = left.merge(
        right,
        on=["regime", "replicate", "seed", "data_hash"],
        suffixes=("_symmetric", "_markov"),
        validate="one_to_one",
    )
    assert len(paired) == 80 and len(traced) == 8
    return native, paired, pd.DataFrame(traced)


def main():
    hashes = inputs_hashes()
    old = json.loads((ALTERNATIVES / "evidence.json").read_text())
    source = json.loads((ALTERNATIVES / "source_change_review.json").read_text())
    assert source["archive_sha256"] == digest(ALTERNATIVES / "source_during_run.py.txt")
    assert (
        source["archive_sha256"] == old["hashes"][str((ALTERNATIVES / "run.py").relative_to(ROOT))]
    )
    assert source["current_sha256"] == digest(ALTERNATIVES / "run.py")
    enriched, recovery, spectra = score_magnitudes()
    native, paired, traces = native_magnitudes(spectra)
    tables = {
        "alternatives_with_magnitudes.csv.gz": enriched,
        "factor_pairs.csv": recovery,
        "native_rows.csv": native,
        "native_pairs.csv": paired,
        "root_trace_scales.csv": traces,
    }
    for name, table in tables.items():
        table.to_csv(HERE / name, index=False)
    assert hashes == inputs_hashes()
    (HERE / "manifest.json").write_text(
        json.dumps(
            {
                "regimes": REGIMES,
                "datasets": len(spectra),
                "inputs_sha256": hashes,
                "outputs": {
                    name: {"rows": len(table), "sha256": digest(HERE / name)}
                    for name, table in tables.items()
                },
                "clustering_refits": 0,
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
