#!/usr/bin/env python3
"""Frozen isolated Wald-rank comparison; no fitted trees or empirical calibration.

Run from the repository root with its existing .venv interpreter. This script
compares the committed prefix rule with the current complete-eigenvalue-group
rule using matched isotropic Gaussian contrasts. All writes stay beside it.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import platform
import subprocess
import sys
import types
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta, binom, chi2

OUTPUT = Path(__file__).resolve().parent
ROOT = OUTPUT.parents[2]
sys.path.insert(0, str(ROOT))

from benchmarks.diagnostics.calibration.reporting import write_diagnostic_bundle
from tree_break_selection.eigenvalue_groups import extend_rank_through_ties
from tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald.projected_wald_kernel import (
    run_projected_wald_kernel,
)

N = 1_000_000
ALPHAS = (0.001, 0.005, 0.01, 0.05)
SEEDS = (1907, 2026092801, 2026092802)
LEGACY_REF = "f34f59cacd09b19c34fafa3a4a6e783dd3973f72"
VERSIONS = ("legacy_prefix", "current_complete_groups")
SELECTIONS = ("fixed_cutoff", "fixed_full", "adaptive_90")
PACKAGE = "tree_break_selection.hierarchy_analysis.statistics.projection.projected_wald"
KERNEL_PATH = PACKAGE.replace(".", "/") + "/projected_wald_kernel.py"
BASIS_PATH = PACKAGE.replace(".", "/") + "/projected_wald_projection_basis.py"
SOURCE_PATHS = (
    KERNEL_PATH,
    BASIS_PATH,
    PACKAGE.replace(".", "/") + "/projected_wald_reference_distribution.py",
    "tree_break_selection/eigenvalue_groups.py",
)


@dataclass(frozen=True)
class Spectrum:
    name: str
    eigenvalues: tuple[float, ...]
    fixed_cutoff: int
    rotation_start: int | None

    @property
    def dimension(self) -> int:
        return len(self.eigenvalues)


SPECTRA = (
    Spectrum("fully_tied_2", (1.0, 1.0), 1, 0),
    Spectrum("distinct_2", (2.0, 1.0), 1, None),
    Spectrum("partial_tie_4", (4.0, 2.0, 2.0, 1.0), 2, 1),
)


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git(*args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True)


def load_legacy_kernel():
    """Use the actual recorded baseline kernel and its own baseline basis helper."""
    sources = {path: git("show", f"{LEGACY_REF}:{path}") for path in (KERNEL_PATH, BASIS_PATH)}
    modules = {}
    for label, path in (("basis", BASIS_PATH), ("kernel", KERNEL_PATH)):
        name = f"{PACKAGE}._wald_rank_audit_legacy_{label}"
        module = types.ModuleType(name)
        module.__package__ = PACKAGE
        sys.modules[name] = module
        exec(compile(sources[path], f"git:{LEGACY_REF}:{path}", "exec"), module.__dict__)
        modules[label] = module
    modules["kernel"].build_pca_projection_basis = modules["basis"].build_pca_projection_basis
    return modules["kernel"].run_projected_wald_kernel, sources


def nominal_rank(spectrum: Spectrum, selection: str) -> int:
    return spectrum.dimension if selection == "fixed_full" else spectrum.fixed_cutoff


def batch(z, spectrum, selection, version, basis):
    """Vectorize the inspected kernel, including its low-energy fallback."""
    projected = z @ basis.T
    energy = projected * projected
    cumulative = np.cumsum(energy, axis=1)
    if selection == "adaptive_90":
        total = np.sum(energy, axis=1)
        rank = np.sum(cumulative < (0.90 * total)[:, None], axis=1) + 1
        rank[total <= 1e-24] = 1
        rank = np.minimum(rank, spectrum.dimension)
    else:
        rank = np.full(len(z), nominal_rank(spectrum, selection), dtype=np.int64)
    if version == "current_complete_groups":
        lookup = np.array(
            [
                extend_rank_through_ties(np.array(spectrum.eigenvalues), k)
                for k in range(spectrum.dimension + 1)
            ]
        )
        rank = lookup[rank]
    statistic = cumulative[np.arange(len(z)), rank - 1]
    return rank.astype(np.uint8), statistic, chi2.sf(statistic, rank)


def rotation(spectrum):
    basis = np.eye(spectrum.dimension)
    if spectrum.rotation_start is None:
        basis[0, 0] = -1.0  # Sign symmetry only: do not mix distinct eigenvalues.
    else:
        start = spectrum.rotation_start
        angle = np.pi / 4
        basis[start : start + 2, start : start + 2] = [
            [np.cos(angle), np.sin(angle)],
            [-np.sin(angle), np.cos(angle)],
        ]
    return basis


def exact_interval(count, n, level=0.95):
    tail = (1.0 - level) / 2.0
    return (
        0.0 if count == 0 else float(beta.ppf(tail, count, n - count + 1)),
        1.0 if count == n else float(beta.ppf(1 - tail, count + 1, n - count)),
    )


def boundary_vectors(dimension):
    samples = [("zero", np.zeros(dimension))]
    for axis in range(dimension):
        z = np.zeros(dimension)
        z[axis] = 3.0
        samples.append((f"axis_{axis}", z))
    for size in (1e-14, 1e-12, 2e-12):
        samples.append((f"small_{size}", np.full(dimension, size)))
    for offset in (-1e-10, 0.0, 1e-10):
        z = np.zeros(dimension)
        z[0], z[-1] = 3.0 + offset, 1.0
        samples.append((f"first_energy_boundary_{offset}", z))
        if dimension == 4:
            z = np.array([0.0, 3.0 + offset, 0.0, 1.0])
            samples.append((f"tie_group_energy_boundary_{offset}", z))
    return samples


def check_kernel(z, spectrum, selection, version, basis, basis_name, kernel):
    samples = [(f"representative_{i}", z[i]) for i in range(64)]
    samples += boundary_vectors(spectrum.dimension)
    values = np.array([value for _label, value in samples])
    ranks, stats, ps = batch(values, spectrum, selection, version, basis)
    rows = []
    for index, (label, value) in enumerate(samples):
        result = kernel(
            value,
            spectral_k=nominal_rank(spectrum, selection),
            pca_projection=basis,
            pca_eigenvalues=np.array(spectrum.eigenvalues),
            adaptive_dimension_energy_fraction=0.9 if selection == "adaptive_90" else None,
        )
        matches = bool(
            int(ranks[index]) == result.projection_dimension
            and np.isclose(stats[index], result.statistic, atol=2e-13, rtol=2e-13)
            and np.isclose(ps[index], result.p_value, atol=2e-14, rtol=2e-13)
        )
        rows.append(
            {
                "spectrum": spectrum.name,
                "selection": selection,
                "version": version,
                "basis": basis_name,
                "sample": label,
                "contrast": json.dumps(value.tolist()),
                "batch_rank": int(ranks[index]),
                "kernel_rank": result.projection_dimension,
                "batch_statistic": float(stats[index]),
                "kernel_statistic": result.statistic,
                "batch_p_value": float(ps[index]),
                "kernel_p_value": result.p_value,
                "matches": matches,
            }
        )
    return rows


def main():
    start = datetime.now(UTC).isoformat()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    initial_sources = {path: (ROOT / path).read_bytes() for path in SOURCE_PATHS}
    legacy_kernel, legacy_sources = load_legacy_kernel()
    design = {
        "schema_version": "isolated_wald_rank_calibration/v1",
        "created_at": start,
        "contrast_count_per_seed_dimension": N,
        "seeds": SEEDS,
        "dimension_4_seed_offset": 4_000_000_000,
        "seed_roles": {
            "1907": "historical replay",
            "2026092801": "held-out",
            "2026092802": "held-out",
        },
        "null": "Independent N(0,I_d) contrasts; eigenvalues label selection groups only.",
        "spectra": [vars(s) for s in SPECTRA],
        "selections": SELECTIONS,
        "versions": VERSIONS,
        "alphas": ALPHAS,
        "legacy_ref": LEGACY_REF,
        "current_source_hashes": {
            path: digest(content) for path, content in initial_sources.items()
        },
        "legacy_source_hashes": {
            path: digest(content.encode()) for path, content in legacy_sources.items()
        },
        "git_head": git("rev-parse", "HEAD").strip(),
        "git_status_at_start": git("status", "--short").splitlines(),
        "script_sha256": digest(Path(__file__).read_bytes()),
        "python": sys.version,
        "platform": platform.platform(),
        "dependencies": {
            name: importlib.metadata.version(name) for name in ("numpy", "scipy", "pandas")
        },
        "command": ".venv/bin/python reports/eigenspace_reassessment_20260928/wald_rank/run_calibration.py",
        "inference": "Exact binomial upper-tail tests of rejection probability <= alpha; Holm FWER .05 over all216 seed/spectrum/selection/version/alpha cells. Pointwise95 and Bonferroni simultaneous95 Clopper-Pearson intervals. Held-out current-rule cells are the primary reassessment; replay is descriptive.",
        "frozen_before_results": True,
        "limitations": [
            "Isolated selected-rank null experiment, not a full pipeline calibration study.",
            "No covariance/PCA fitting, tree construction, diffusion, branch lengths, edge selection, empirical inflation, sibling BH or traversal.",
            "Shared contrasts across spectra and rules are paired comparisons; multiplicity correction does not require cell independence.",
            "Fixed spectra are metadata, not covariance eigenvalues of the isotropic null.",
            "No rule or parameter is tuned to held-out outcomes.",
        ],
    }
    (OUTPUT / "design.json").write_text(json.dumps(design, indent=2) + "\n")
    snapshot = OUTPUT / "source_snapshot"
    snapshot.mkdir(exist_ok=True)
    for path, source in initial_sources.items():
        (snapshot / ("current_" + Path(path).name + ".txt")).write_bytes(source)
    for path, source in legacy_sources.items():
        (snapshot / ("legacy_" + Path(path).name + ".txt")).write_text(source)

    count_rows, rank_rows, paired_rows, rotation_rows, kernel_rows = [], [], [], [], []
    outcome_index = []
    expected_cells = len(SEEDS) * len(SPECTRA) * len(SELECTIONS) * len(VERSIONS) * len(ALPHAS)
    for seed in SEEDS:
        draws = {
            d: np.random.default_rng(seed if d == 2 else seed + 4_000_000_000).normal(size=(N, d))
            for d in (2, 4)
        }
        packed = {}
        for spectrum in SPECTRA:
            z = draws[spectrum.dimension]
            for selection in SELECTIONS:
                outcomes = {}
                for version in VERSIONS:
                    rank, statistic, p_values = batch(
                        z, spectrum, selection, version, np.eye(spectrum.dimension)
                    )
                    outcomes[version] = (rank, p_values)
                    name = f"{spectrum.name}__{selection}__{version}"
                    packed[name + "__rank"] = rank
                    counts = []
                    for alpha in ALPHAS:
                        rejected = p_values <= alpha
                        count = int(rejected.sum())
                        counts.append(count)
                        packed[name + "__alpha_" + str(alpha)] = np.packbits(
                            rejected, bitorder="little"
                        )
                        ci = exact_interval(count, N)
                        simultaneous = exact_interval(count, N, 1 - 0.05 / expected_cells)
                        count_rows.append(
                            {
                                "seed": seed,
                                "role": "replay" if seed == 1907 else "held_out",
                                "spectrum": spectrum.name,
                                "dimension": spectrum.dimension,
                                "selection": selection,
                                "nominal_rank": nominal_rank(spectrum, selection),
                                "version": version,
                                "alpha": alpha,
                                "n_contrasts": N,
                                "rejections": count,
                                "rejection_rate": count / N,
                                "rate_minus_alpha": count / N - alpha,
                                "cp95_low": ci[0],
                                "cp95_high": ci[1],
                                "simultaneous_cp95_low": simultaneous[0],
                                "simultaneous_cp95_high": simultaneous[1],
                                "raw_upper_binomial_p": float(binom.sf(count - 1, N, alpha)),
                            }
                        )
                    for value, count in zip(*np.unique(rank, return_counts=True), strict=True):
                        rank_rows.append(
                            {
                                "seed": seed,
                                "spectrum": spectrum.name,
                                "selection": selection,
                                "version": version,
                                "rank": int(value),
                                "count": int(count),
                            }
                        )
                    outcome_index.append(
                        {
                            "archive": f"replicate_outcomes_seed_{seed}.npz",
                            "key_prefix": name,
                            "seed": seed,
                            "n_contrasts": N,
                            "dimension": spectrum.dimension,
                        }
                    )
                    if seed == 1907:
                        kernel = (
                            legacy_kernel
                            if version == "legacy_prefix"
                            else run_projected_wald_kernel
                        )
                        for basis_name, basis in (
                            ("identity", np.eye(spectrum.dimension)),
                            ("symmetry", rotation(spectrum)),
                        ):
                            kernel_rows.extend(
                                check_kernel(
                                    z, spectrum, selection, version, basis, basis_name, kernel
                                )
                            )
                        r0, t0, p0 = batch(
                            z[:100000], spectrum, selection, version, np.eye(spectrum.dimension)
                        )
                        r1, t1, p1 = batch(
                            z[:100000], spectrum, selection, version, rotation(spectrum)
                        )
                        row = {
                            "spectrum": spectrum.name,
                            "selection": selection,
                            "version": version,
                            "n_contrasts": 100000,
                            "rotation": "45 degrees within repeated eigenspace"
                            if spectrum.rotation_start is not None
                            else "first eigenvector sign flip",
                            "rank_mismatches": int(np.sum(r0 != r1)),
                            "max_abs_statistic_difference": float(np.max(np.abs(t0 - t1))),
                            "max_abs_p_difference": float(np.max(np.abs(p0 - p1))),
                        }
                        for alpha in ALPHAS:
                            row["decision_mismatches_alpha_" + str(alpha)] = int(
                                np.sum((p0 <= alpha) != (p1 <= alpha))
                            )
                        rotation_rows.append(row)
                    print(seed, spectrum.name, selection, version, "rejections", counts, flush=True)
                old_r, old_p = outcomes["legacy_prefix"]
                new_r, new_p = outcomes["current_complete_groups"]
                for alpha in ALPHAS:
                    old, new = old_p <= alpha, new_p <= alpha
                    paired_rows.append(
                        {
                            "seed": seed,
                            "spectrum": spectrum.name,
                            "selection": selection,
                            "alpha": alpha,
                            "n_contrasts": N,
                            "old_only_rejected": int(np.sum(old & ~new)),
                            "new_only_rejected": int(np.sum(new & ~old)),
                            "both_rejected": int(np.sum(old & new)),
                            "neither_rejected": int(np.sum(~old & ~new)),
                            "rank_changed": int(np.sum(old_r != new_r)),
                        }
                    )
        np.savez_compressed(OUTPUT / f"replicate_outcomes_seed_{seed}.npz", **packed)
    counts = pd.DataFrame(count_rows)
    raw = counts["raw_upper_binomial_p"].to_numpy()
    order = np.argsort(raw)
    corrected = np.minimum(
        1.0, np.maximum.accumulate(raw[order] * (len(raw) - np.arange(len(raw))))
    )
    holm = np.empty_like(corrected)
    holm[order] = corrected
    counts["holm_adjusted_p"] = holm
    counts["holm_excess_rejection_flag"] = holm <= 0.05
    assert len(counts) == expected_cells == 216
    kernel_table = pd.DataFrame(kernel_rows)
    mismatches = kernel_table.loc[~kernel_table["matches"]]
    tables = {
        "counts": counts,
        "rank_counts": pd.DataFrame(rank_rows),
        "paired_decisions": pd.DataFrame(paired_rows),
        "rotation_audit": pd.DataFrame(rotation_rows),
        "kernel_checks": kernel_table,
        "outcome_index": pd.DataFrame(outcome_index),
    }
    current_unchanged = all(
        (ROOT / path).read_bytes() == source for path, source in initial_sources.items()
    )
    bundle = write_diagnostic_bundle(
        output_dir=OUTPUT,
        tables=tables,
        filenames={name: name + ".csv" for name in tables},
        manifest={
            **design,
            "finished_at": datetime.now(UTC).isoformat(),
            "all_kernel_checks_match": bool(mismatches.empty),
            "kernel_check_count": len(kernel_table),
            "production_sources_unchanged_during_run": current_unchanged,
            "multiplicity_family_size": len(counts),
            "raw_outcome_format": "npz rank uint8 and bitpacked rejection flags; replicate identifier is zero-based vector index; np.unpackbits(...,bitorder='little')[:n_contrasts]",
            "outcome_archive_sha256": {
                f"replicate_outcomes_seed_{seed}.npz": digest(
                    (OUTPUT / f"replicate_outcomes_seed_{seed}.npz").read_bytes()
                )
                for seed in SEEDS
            },
        },
    )
    if not mismatches.empty:
        print(mismatches.to_string(index=False), flush=True)
        raise RuntimeError("Vectorized/kernel parity failed; inspect saved kernel_checks.csv")
    if not current_unchanged:
        raise RuntimeError("Production sources changed during run")
    print(json.dumps({name: str(path) for name, path in bundle.items()}, indent=2))


if __name__ == "__main__":
    main()
