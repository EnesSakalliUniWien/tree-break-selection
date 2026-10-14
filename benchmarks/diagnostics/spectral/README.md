# Spectral Diagnostics

This directory contains diagnostics for projection-dimension and
Marchenko-Pastur threshold behavior.

- `mp/` owns MP-law checks, dimension sweeps, and backend profiling.
- `adaptive_cosine/` owns KAK, lens, and adaptive-cosine probes.
- `stability/` owns covariance-axis and tree-strategy stability panels.
- `vandermonde.py` provides guarded PSD Toeplitz factorization for diagnostics;
  `apply_vandermonde_diagnostic.py` checks its applicability to the known-branch
  covariance experiment and records full-rank nonuniqueness.

Entrypoint:

```bash
uv run python -m benchmarks.diagnostics.spectral.mp.compare_mp_dimension_contracts
```

Backend profiling:

```bash
uv run python -m benchmarks.diagnostics.spectral.mp.profile_spectral_backends
```

Local identity-MP law screen:

```bash
uv run python -m benchmarks.diagnostics.spectral.mp.local_mp_identity_law_diagnostic
```

The local identity-MP screen is descriptive only. It writes case-level MP
departure summaries, node-level spectra, selected-tree spectral-law covariate
relationships, and a separate categorical extreme-node table. It must not be
used as a production threshold or fallback calibration path.

## Eigenvalue block comparisons

`adaptive_cosine/compare_spectral_blocks.py` compares the existing adaptive
cosine blocks on matched benchmark data. It records eigenvalue gaps, energy,
sample-distance agreement, rooted tree differences, TBS assignments and
explicit unavailable outcomes across the seven linkage methods and neighbor
joining, using raw coordinates and adaptive diffusion separately.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MPLBACKEND=Agg uv run --no-sync python -m benchmarks.diagnostics.spectral.adaptive_cosine.compare_spectral_blocks \
  --suite method_proof --output-dir reports/new-spectral-block-comparison
```

The output directory must be new. Defaults retain a target of 40 positive
modes, preserve numerical ties, and use the standard TBS runner with edge
alpha 0.001, sibling alpha 0.01 and internal support thresholds unenforced.
`--enforce-internal-support-thresholds` enables the strict support check.
Ground-truth labels are used only for evaluation. NNLS branch-time variants
are outside this distance/topology comparison.
The analysis is descriptive and does not establish calibration or select a
winning block or method.
