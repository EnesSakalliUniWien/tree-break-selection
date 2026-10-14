# Repeated-eigenspace rank checks and repair — 2026-09-28

The checked diffusion, local PCA and Wald rank boundaries now retain complete
numerical eigenvalue groups. The reproduced sample/feature permutations and
rotations within repeated eigenspaces no longer change their geometry or Wald
decisions. This repairs rank-selection geometry; adaptive-rank and selected-tree
calibration remain unvalidated.

## Before and after

The same review probes were rerun against the repaired source. Exact results,
source hashes and verification counts are in
[rank_invariance_repair.json](/Users/berksakalli/Projects/kl-te-cluster/reports/mathematical_construction_review_20260928/rank_invariance_repair.json).
The original [findings.json](/Users/berksakalli/Projects/kl-te-cluster/reports/mathematical_construction_review_20260928/findings.json)
is preserved as the snapshot before this repair.

| Check | Before | After |
|---|---|---|
| Six-cycle diffusion, target rank one, sample permutation | Distance change 0.21650635094611065 | Retains both tied directions; change 9.71445146547012e-16 |
| Covariance spectrum 0.5, 0.5, target rank one, feature permutation | Projector difference sqrt(2) | Retains rank two; projector difference zero |
| Adaptive Wald, contrast (3,0), basis rotation 45 degrees | p changes 0.002699796 to 0.011108997 | Both p-values 0.011108997, rank two, statistic nine |
| Dense continuous 20-by-2 whitening example | MP raw count zero, test rank one | MP raw count still zero, test rank two, complete repeated space retained |

## Implementation contract

- [eigenvalue_groups.py](/Users/berksakalli/Projects/kl-te-cluster/tree_break_selection/eigenvalue_groups.py)
  extracts the existing adaptive-cosine grouping rule for shared use: adjacent
  eigenvalues join transitively under relative tolerance `1e-10` and absolute
  tolerance `1e-12 * spectral_scale`. These are numerical groups, not exact
  multiplicity certificates.
- Cosine, diffusion, PCA and Wald rank arguments are targets that extend to the
  end of the boundary group. PCA/MP expansion remains bounded by the existing
  numerically positive available rank; zero target/zero signal with floor zero
  remains zero.
- Diffusion removes the first mode only when its eigenvalue is isolated. A
  repeated leading group stays whole, including on disconnected graphs. This
  follows the existing cosine rule for avoiding arbitrary common-mode isolation.
  The symmetric metric and negative-eigenvalue clipping are unchanged.
- MP metadata reports the extended **test dimension**, while the **raw MP
  signal count** retains its original threshold meaning. Wald degrees of freedom
  follow the actual returned basis dimension.
- PCA returns its actual positive eigenvalues. The former `1e-12` output floor
  could otherwise turn distinct small eigenvalues into artificial ties for the
  downstream Wald helper. The separate numerical protections in dual recovery
  are retained.
- A Wald parent basis that is known to omit part of the required tied group is
  rejected explicitly; the helper cannot invent missing eigenvectors.
- The application CLI help and the space-separation/statistics documentation
  describe rank extension. Diffusion metadata already records the actual rank.

## Verification

The 31 existing targeted baseline cases passed before editing. The first 21 new
regression cases produced **18 failures and three passing controls**. The final
suite has **25 passing regression cases** in
[55_test_eigenspace_invariance.py](/Users/berksakalli/Projects/kl-te-cluster/tests/statistics/55_test_eigenspace_invariance.py).

Coverage includes exact and numerical ties, transitive groups, sample and feature
permutations, rotations of tied bases, primal and dual PCA, repeated leading
diffusion groups, distinct-spectrum controls, positive-rank limits, small positive
eigenvalues, MP threshold/floor metadata, incomplete parent bases, and fixed and
adaptive Wald decisions. Two older Wald tests now use distinct eigenvalues to
continue testing their intended nominal-rank behavior; tied cases are explicit
in the new suite.

```bash
.venv/bin/python -m pytest -q -p no:cacheprovider tests/statistics/55_test_eigenspace_invariance.py
make check
```

Full `make check` passed **693 tests** across its seven stages, with stage counts
192, 299, 1, 14, 144, 40 and 3. Ruff, dependency audit, Vulture and wiki lint also
passed. A separate read-only mathematical/integration review found no material
defect in the changes. The reviewer inspected the source/tests and did not run
an additional test suite.

## Remaining limitations

Retaining complete groups can increase dimension, memory use and runtime,
particularly for flat spectra after dense continuous self-whitening. The checks
establish the documented examples and integration behavior, not a bound on every
dataset's numerical sensitivity or resource use.

No full benchmark panel or new full-pipeline null study was run. The old
one-million-contrast result belongs to the pre-repair source snapshot; it must
not be presented as a post-repair error estimate. Selecting rank from the tested
contrast still needs a justified reference law when distinct eigenvalue groups
allow that rank to vary. The whitening/MP interpretation also remains an open
question.

The separate unequal-rank principal-angle diagnostic, Hamming self-neighbor
issue and broader documentation contradictions were not part of this rank
repair. Existing worktree changes and historical numerical evidence were
preserved. No commit or push was made.
