"""Known-identity Gaussian selected-root bounds; a diagnostic, not production traversal.

Run with ``uv run --no-sync python -m
benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_diagnostic
--output-dir PATH``. Unresolved cases remain in the denominator of the report.
"""

import argparse
import hashlib
import json
import sys
from decimal import Decimal
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd

from benchmarks.diagnostics.calibration.reporting import (
    print_diagnostic_output_paths,
    write_diagnostic_bundle,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_cells import (
    analyze_diffusion_root_cells,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_mass import (
    bounded_root_gate,
    chi_interval_mass,
)
from benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_replay import (
    make_diffusion_root_path,
)


def diagnose_diffusion_root(
    values,
    *,
    alpha=0.01,
    anchor_ratios=(1.0,),
    subdivision_depth=0,
    observed_half_width=None,
    expansion_steps=0,
):
    """Replay, certify and gate one root without any fitted test nuisance parameters."""
    threshold = Decimal(str(alpha))
    if not threshold.is_finite() or not 0 < threshold <= 1:
        raise ValueError("alpha must be finite in (0,1]")
    shape = np.shape(values)
    if len(shape) != 2:
        raise ValueError("values must be an n-by-p matrix")
    root = dict(
        sample_size=shape[0],
        feature_dimension=shape[1],
        covariance="known_identity",
        projection="full_feature",
        branch_length=0.0,
        inflation=1.0,
        contrast_energy_rank_selection=False,
        alpha=str(threshold),
        passthrough=False,
        input_data_sha256=hashlib.sha256(np.asarray(values, dtype=float).tobytes()).hexdigest(),
        observed_half_width=observed_half_width,
        expansion_steps=expansion_steps,
    )
    diagnostic_hashes = {
        str(source.resolve()): hashlib.sha256(source.read_bytes()).hexdigest()
        for source in Path(__file__).parent.glob("diffusion_root_*.py")
    }
    try:
        path = make_diffusion_root_path(values, sigma=1.0)
    except ValueError as error:
        root.update(
            gate_status="unsupported",
            clusters=1,
            p_lower="0",
            p_upper="1",
            selection_mass_lower="0",
            selection_mass_upper="1",
            error=str(error),
        )
        return dict(
            root=root,
            cells=[],
            replay=[],
            source_hashes=diagnostic_hashes,
            arithmetic_contract=None,
        )
    analysis = analyze_diffusion_root_cells(
        path,
        anchor_ratios=anchor_ratios,
        subdivision_depth=subdivision_depth,
        observed_half_width=observed_half_width,
        expansion_steps=expansion_steps,
    )
    bounds, gate = bounded_root_gate(analysis, path, alpha)
    root.update(
        canonical_data_sha256=hashlib.sha256(path.values.tobytes()).hexdigest(),
        observed_radius=path.observed_radius,
        observed_radius_bounds=json.dumps([str(v) for v in analysis.observed_radius_bounds]),
        observed_history_certified=analysis.selection_certified,
        gate_status=gate.status,
        clusters=gate.clusters,
        p_lower=str(bounds.lower),
        p_upper=str(bounds.upper),
        selection_mass_lower=str(bounds.selection_mass_lower),
        selection_mass_upper=str(bounds.selection_mass_upper),
        unresolved_below_lower=str(bounds.unresolved_below[0]),
        unresolved_below_upper=str(bounds.unresolved_below[1]),
        unresolved_above_lower=str(bounds.unresolved_above[0]),
        unresolved_above_upper=str(bounds.unresolved_above[1]),
    )
    cells = []
    for cell in analysis.cells:
        mass = chi_interval_mass(cell.lo, cell.hi, path.dimension)
        cells.append(
            dict(
                lo=str(cell.lo),
                hi=str(cell.hi),
                status=cell.status,
                proof=cell.proof,
                mass_lower=str(mass[0]),
                mass_upper=str(mass[1]),
                details=json.dumps(cell.details, sort_keys=True, default=str),
            )
        )
    replay = []
    for ratio in anchor_ratios:
        point = path.evaluate(float(ratio) * path.observed_radius)
        replay.append(
            dict(
                radius_ratio=float(ratio),
                radius=point.radius,
                status=point.status,
                minimum_point_margin=float(point.margins.min()) if len(point.margins) else None,
                signature=json.dumps(point.signature),
                error=point.error,
                kde_dimension=point.state.get("kde_dimension"),
                kde_bgh_index=point.state.get("kde_bgh_index"),
                outer_epsilon=point.state.get("outer_epsilon"),
                retained_dimension=point.state.get("retained_dimension"),
                positive_self_scaled_distance_count=point.state.get(
                    "positive_self_scaled_distance_count"
                ),
                actual_positive_scaled_sq_median=point.state.get(
                    "actual_positive_scaled_sq_median"
                ),
                off_diagonal_positive_scaled_sq_median=point.state.get(
                    "off_diagonal_positive_scaled_sq_median"
                ),
                interval_certificate=False,
            )
        )
    hashes = dict(analysis.source_hashes)
    for source in (Path(__file__), Path(chi_interval_mass.__code__.co_filename)):
        hashes[str(source.resolve())] = hashlib.sha256(source.read_bytes()).hexdigest()
    return dict(
        root=root,
        cells=cells,
        replay=replay,
        source_hashes=hashes,
        arithmetic_contract=analysis.arithmetic_contract,
    )


def write_root_diagnostics(results, output_dir):
    """Write separate point observations, interval certificates and bounded root decisions."""
    if not results:
        raise ValueError("at least one diagnostic result is required")
    tables = {}
    for key in ("root", "cells", "replay"):
        rows = []
        for replicate, result in enumerate(results):
            records = [result[key]] if key == "root" else result[key]
            rows.extend(dict(replicate=replicate, **row) for row in records)
        table = pd.DataFrame(rows)
        if table.empty:
            table = pd.DataFrame(columns=["replicate", "status"])
        tables["roots" if key == "root" else key] = table
    hashes = {}
    for result in results:
        for path, digest in result["source_hashes"].items():
            if path in hashes and hashes[path] != digest:
                raise ValueError("source changed during the diagnostic panel")
            hashes[path] = digest
    contracts = {r["arithmetic_contract"] for r in results if r["arithmetic_contract"]}
    if len(contracts) > 1:
        raise ValueError("inconsistent arithmetic contracts")
    counts = tables["roots"].gate_status.value_counts().to_dict()
    return write_diagnostic_bundle(
        output_dir=Path(output_dir),
        tables=tables,
        filenames={name: f"{name}.csv" for name in tables},
        manifest=dict(
            scope="diagnostic_root_only",
            conditioning="complete_observed_merge_history",
            null="independent Gaussian rows, known identity covariance, eta.T @ E[X] = 0",
            arithmetic_contract=next(iter(contracts), "no supported observed constructor"),
            source_hashes=hashes,
            python=sys.version.split()[0],
            versions={
                name: version(name) for name in ("numpy", "scipy", "scikit-learn", "pydiffmap")
            },
            root_counts=counts,
            limitations=[
                "Point replays are observations, not interval certificates.",
                "The native identity diagonal is zeroed before epsilon; interval certificates still concern exact-real arithmetic.",
                "Unresolved finite intervals and infinite tails retain their probability mass.",
                "Root rejection opens two children only; descendant traversal is not calibrated.",
                "A small diagnostic panel does not establish empirical calibration.",
            ],
        ),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--sample-size", type=int, default=9)
    parser.add_argument("--features", type=int, default=1)
    parser.add_argument("--seed", type=int, default=20261105)
    parser.add_argument("--subdivision-depth", type=int, default=0)
    parser.add_argument("--observed-half-width", type=float)
    parser.add_argument("--expansion-steps", type=int, default=0)
    args = parser.parse_args()
    if min(args.replicates, args.sample_size, args.features) < 1:
        parser.error("replicates, sample-size and features must be positive")
    results = []
    for replicate in range(args.replicates):
        seed = args.seed + replicate
        values = np.random.default_rng(seed).normal(size=(args.sample_size, args.features))
        result = diagnose_diffusion_root(
            values,
            anchor_ratios=(1.0, 1.1, 1.2),
            subdivision_depth=args.subdivision_depth,
            observed_half_width=args.observed_half_width,
            expansion_steps=args.expansion_steps,
        )
        result["root"]["seed"] = seed
        results.append(result)
    print_diagnostic_output_paths(write_root_diagnostics(results, args.output_dir))


if __name__ == "__main__":
    main()
