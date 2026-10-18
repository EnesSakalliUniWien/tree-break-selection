"""Apply guarded Vandermonde decomposition to the known-branch covariance study."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import subspace_angles

from benchmarks.diagnostics.calibration.reporting import write_diagnostic_bundle
from benchmarks.diagnostics.calibration.sibling.known_branch_skew_diagnostic import (
    LEAF_NAMES,
    covariance_matrix,
    population_means,
)
from benchmarks.diagnostics.spectral.vandermonde import (
    decompose_psd_toeplitz,
    sinusoidal_atoms,
    toeplitz_residual,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--protocol", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        raise FileExistsError(args.output_dir)
    root = Path(__file__).resolve().parents[3]
    paths = [
        Path(__file__).resolve(),
        Path(__file__).with_name("vandermonde.py"),
        args.protocol.resolve(),
        args.reference.resolve() / "manifest.json",
        root / "benchmarks/diagnostics/calibration/sibling/known_branch_skew_diagnostic.py",
    ]
    hashes = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
    args.output_dir.mkdir(parents=True)
    (args.output_dir / "source_hashes_before.json").write_text(json.dumps(hashes, indent=2) + "\n")
    manifest = json.loads((args.reference / "manifest.json").read_text())
    config = manifest["config"]
    p, n = config["dimension"], config["n_per_group"]
    reference = pd.read_csv(args.reference / "known_branch_selected_root_rows.csv").set_index(
        ["covariance_profile", "pattern", "replicate"]
    )
    seeds = iter(np.random.SeedSequence(config["base_seed"]).spawn(len(reference)))
    cases = []
    structure_rows = []
    for profile in config["covariance_profiles"]:
        sigma = covariance_matrix(profile, p)
        for anchor in (0.0, 0.137):
            cases.append((f"population_{profile}", sigma, anchor))
        chol = np.linalg.cholesky(sigma)
        for pattern in config["patterns"]:
            means = population_means(pattern, p)
            for replicate in range(config["replicates"]):
                rng = np.random.default_rng(next(seeds))
                data = np.vstack(
                    [rng.multivariate_normal(means[leaf], sigma, size=n) for leaf in LEAF_NAMES]
                )
                digest = hashlib.sha256(np.ascontiguousarray(data).tobytes()).hexdigest()
                assert reference.loc[(profile, pattern, replicate), "dataset_sha256"] == digest
                centered = data - data.mean(axis=0)
                sample = centered.T @ centered / (len(data) - 1)
                known_white = np.linalg.solve(chol, centered.T).T
                own_white = np.linalg.solve(np.linalg.cholesky(sample), centered.T).T
                matrices = {
                    "sample_pooled": sample,
                    "oracle_whitened_PCA": known_white.T @ known_white / len(data),
                    "sample_whitened_PCA": own_white.T @ own_white / len(data),
                }
                for name, matrix in matrices.items():
                    residual = toeplitz_residual(matrix)
                    values = np.linalg.eigvalsh(matrix)
                    structure_rows.append(
                        {
                            "profile": profile,
                            "pattern": pattern,
                            "replicate": replicate,
                            "dataset_sha256": digest,
                            "matrix_kind": name,
                            "toeplitz_relative_residual": residual,
                            "positive_definite": bool(values[0] > 0),
                            "applicable_at_tolerance": residual <= 1e-10 and values[0] >= 0,
                        }
                    )
                    if pattern == "global_null" and replicate == 0:
                        cases.append((f"{profile}_{name}_rep0", matrix, 0.0))
    for name, frequencies, weights in (
        ("known_complex_rank2", [0.12, 0.37], [0.7, 1.3]),
        ("known_real_conjugate_rank2", [0.125, 0.875], [0.5, 0.5]),
    ):
        atoms = sinusoidal_atoms(p, np.array(frequencies))
        cases.append((name, (atoms * weights) @ atoms.conj().T, 0.0))
    rows, atoms_rows = [], []
    for name, matrix, anchor in cases:
        row = {
            "case": name,
            "anchor": anchor,
            "toeplitz_relative_residual": toeplitz_residual(matrix),
        }
        try:
            result = decompose_psd_toeplitz(matrix, anchor_frequency=anchor)
            values, vectors = np.linalg.eigh(matrix)
            eigenbasis = vectors[:, -result.numerical_rank :]
            q, _ = np.linalg.qr(result.atoms)
            gram = result.atoms.conj().T @ result.atoms / len(matrix)
            reconstruction = (result.atoms * result.weights) @ result.atoms.conj().T
            np.testing.assert_allclose(reconstruction, matrix, atol=1e-9, rtol=1e-9)
            row.update(
                status="ok",
                failure_reason="",
                numerical_rank=result.numerical_rank,
                unique_at_numerical_rank=result.unique_at_numerical_rank,
                relative_reconstruction_error=result.relative_reconstruction_error,
                atom_condition=result.atom_condition,
                max_subspace_angle_degrees=float(np.degrees(subspace_angles(q, eigenbasis)).max()),
                projector_frobenius_error=float(
                    np.linalg.norm(q @ q.conj().T - eigenbasis @ eigenbasis.conj().T)
                ),
                max_atom_normalized_overlap=float(
                    np.max(np.abs(gram - np.eye(result.numerical_rank)))
                ),
            )
            for frequency, weight in zip(result.frequencies, result.weights):
                atoms_rows.append(
                    {
                        "case": name,
                        "anchor": anchor,
                        "frequency": frequency,
                        "positive_weight": weight,
                    }
                )
        except ValueError as exc:
            row.update(status="unsupported", failure_reason=str(exc))
        rows.append(row)
    for path, digest in hashes.items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest, path
    structure = pd.DataFrame(structure_rows)
    summary = structure.groupby(["profile", "pattern", "matrix_kind"], as_index=False).agg(
        n=("applicable_at_tolerance", "size"),
        n_applicable=("applicable_at_tolerance", "sum"),
        min_toeplitz_residual=("toeplitz_relative_residual", "min"),
        max_toeplitz_residual=("toeplitz_relative_residual", "max"),
    )
    outputs = write_diagnostic_bundle(
        output_dir=args.output_dir,
        tables={
            "decompositions": pd.DataFrame(rows),
            "atoms": pd.DataFrame(atoms_rows),
            "structure": structure,
            "structure_summary": summary,
        },
        filenames={
            "decompositions": "decompositions.csv",
            "atoms": "atoms.csv",
            "structure": "structure.csv",
            "structure_summary": "structure_summary.csv",
        },
        manifest={
            "schema_version": "vandermonde_diagnostic/v1",
            "source_hashes_before_and_after_equal": True,
            "source_sha256": hashes,
            "reference": str(args.reference),
            "config": config,
            "toeplitz_relative_tolerance": 1e-10,
            "rank_relative_tolerance": 1e-10,
            "decomposition_source": "https://arxiv.org/abs/1605.02431",
            "interpretation": "Exact PSD Toeplitz factorization only; no fitting or production changes. Full-rank frequencies are nonunique. Numerical rank is not exact certification.",
        },
    )
    print(json.dumps({key: str(value) for key, value in outputs.items()}, indent=2))


if __name__ == "__main__":
    main()
