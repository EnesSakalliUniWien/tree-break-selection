# Isolated projected-Wald rank reassessment

Retaining complete repeated eigenspaces repaired rotation invariance. It also
removed the adaptive-rank excess rejection in the fully tied two-dimensional
example, where the repaired rule always retains both directions. It did **not**
calibrate contrast-selected rank across distinct or partially tied groups.

This is an isolated `N(0,I_d)` contrast experiment. The supplied PCA eigenvalues
identify selection groups; they are not covariance eigenvalues of the simulated
null. There is no estimated covariance, tree, diffusion, branch length, edge
selection, empirical inflation, sibling BH or traversal in this experiment.

## Matched replay

The same one million seed-1907 two-dimensional contrasts used in the earlier
review are reused exactly. Four-dimensional contrasts use seed 4000001907.
Results below use adaptive energy fraction 0.90 at alpha 0.01.

| Spectrum | Legacy prefix rejections | Complete-group rejections | Denominator |
|---|---:|---:|---:|
| `(1,1)` | 14,720 | 10,070 | 1,000,000 |
| `(2,1)` | 14,720 | 14,720 | 1,000,000 |
| `(4,2,2,1)` | 17,575 | 15,434 | 1,000,000 |

The first comparison resolves the original tied-eigenspace counterexample.
The unchanged distinct-spectrum result separates that repair from the remaining
contrast-dependent rank-selection problem.

## Frozen held-out evaluation

Seeds 2026092801 and 2026092802 were specified before their outcomes were read.
Each seed supplies one million contrasts per dimension; dimension four uses
the seed plus 4000000000. Spectra and rules within a dimension use matched
contrasts. No parameter or rule was tuned to these outcomes.

Current adaptive-0.90 results at alpha 0.01, with pointwise exact binomial
95% Clopper-Pearson confidence intervals:

| Spectrum | Seed | Rejection rate | 95% interval |
|---|---:|---:|---:|
| `(1,1)` | 2026092801 | 0.009949 | [0.009755, 0.010145] |
| `(1,1)` | 2026092802 | 0.009993 | [0.009799, 0.010190] |
| `(2,1)` | 2026092801 | 0.014464 | [0.014231, 0.014700] |
| `(2,1)` | 2026092802 | 0.014595 | [0.014361, 0.014832] |
| `(4,2,2,1)` | 2026092801 | 0.015531 | [0.015290, 0.015775] |
| `(4,2,2,1)` | 2026092802 | 0.015342 | [0.015102, 0.015585] |

The complete design evaluates alpha 0.001, 0.005, 0.01 and 0.05; three spectra;
fixed cutoff, fixed full rank and adaptive 0.90; both legacy and current rules;
and three seeds. The fixed cutoff is one for dimension two and two for dimension
four, so the latter also tests a cutoff inside the middle repeated eigenspace.

`counts.csv` contains all 216 cells, pointwise confidence intervals, Bonferroni
simultaneous 95% confidence intervals, and one-sided exact binomial tests of
rejection probability at most alpha. Holm adjustment at family error level 0.05
covers all 216 cells, including dependent matched comparisons and the historical
replay. Exact-tail values stored as zero denote floating-point underflow, not a
literal probability zero; the simultaneous confidence bounds independently
exclude the target in those cells.

All 16 held-out current adaptive cells with distinct or partially tied spectra
show excess rejection after this correction, at every tested alpha and both
seeds. No held-out current fixed-rank control or fully tied adaptive cell is
flagged. Absence of a flag is not a general calibration guarantee.

## Symmetry and numerical verification

On 100,000 matched contrasts per spectrum and rule, rotating by 45 degrees
inside the repeated eigenspace produces no current-rule rank or rejection
changes. The largest current-rule p-value difference is below 2.8e-15.
Distinct eigenvalues are checked only under an eigenvector sign change;
mixing their directions would not be a symmetry of that spectrum.

The legacy adaptive rule changes rank for 40,890 of 100,000 fully tied contrasts
and 14,604 of 100,000 partially tied contrasts under the same rotation. The
current complete-group rule has zero corresponding rank changes. Fixed legacy
cutoffs can change statistics despite retaining the same numerical rank.

All 2,688 vectorized-versus-actual-kernel checks pass. They include ordinary
contrasts, axis vectors, zero and low-energy vectors, and exact and neighboring
energy boundaries, under both the original and rotated bases. The legacy
checks execute the actual kernel and basis helper from commit
`f34f59cacd09b19c34fafa3a4a6e783dd3973f72`. Current source hashes are unchanged
through the run and snapshots are retained.

Before this artifact was created, 33 existing projected-Wald, reference-law and
eigenspace-invariance tests passed. Separately decoding the compressed output
archives reproduces all 216 rejection counts, 108 paired-decision rows and 78
rank-count rows; see `verification.json`.

## Artifacts and reproduction

- `run_calibration.py`: frozen design and reproducible computation.
- `design.json`, `manifest.json`: seeds, versions, source hashes, runtime,
  limitations, inference family and output hashes.
- `counts.csv`: rates, confidence intervals and corrected excess-rejection flags.
- `rank_counts.csv`, `paired_decisions.csv`: rank changes and matched decisions.
- `rotation_audit.csv`, `kernel_checks.csv`: symmetry and implementation checks.
- `replicate_outcomes_seed_*.npz`, `outcome_index.csv`: compact per-contrast
  ranks and rejection outcomes. Contrast identifier is the zero-based vector
  index; decode rejection vectors with
  `np.unpackbits(array, bitorder="little")[:1000000]`. Their generation is fully
  specified by the saved seed, dimension and script.
- `source_snapshot/`: baseline and current source, saved with `.txt` suffixes.
- `verify_outputs.py`, `verification.json`: independent output decoding audit.

From the repository root:

```bash
.venv/bin/python reports/eigenspace_reassessment_20260928/wald_rank/run_calibration.py
.venv/bin/python reports/eigenspace_reassessment_20260928/wald_rank/verify_outputs.py
```

The experiment uses six million generated contrasts across three seeds and two
dimensions, shared across the matched comparisons. It cannot establish or
refute the full pipeline's selected-tail, sibling-FDR or false-partition rate.
