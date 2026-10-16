# Spectral Diagnostics

This directory contains adaptive-cosine, tree-strategy stability, and
Vandermonde diagnostics.

- `adaptive_cosine/` owns KAK, lens, and adaptive-cosine probes.
- `stability/` owns the tree-strategy stability panel.
- `vandermonde.py` provides guarded PSD Toeplitz factorization for diagnostics;
  `apply_vandermonde_diagnostic.py` checks its applicability to the known-branch
  covariance experiment and records full-rank nonuniqueness.

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
