---
title: Cancer PC Surrogate Evidence
type: source
status: reviewed
updated: 2026-11-06
sources:
  - raw/inbox/choir-cell-lines/tbs-merged-rerun-20261106.md
  - raw/inbox/cancer-pc-surrogate-review-20261104.md
  - raw/inbox/cancer-pc-surrogate-experiment-20261104.py.txt
  - raw/assets/cancer-pc-surrogate-20261104-final/manifest.json
  - raw/assets/cancer-pc-surrogate-20261104-final/results.json
  - raw/assets/cancer-pc-surrogate-20261104-final/summary.csv
  - raw/inbox/choir-cell-lines/tbs-root-groups.md
  - raw/inbox/choir-cell-lines/tbs-diagnosis.md
  - benchmarks/shared/generators/cancer_pc_surrogate.py
  - benchmarks/shared/cases/method_proof.py
  - tests/core/test_cancer_pc_surrogate.py
  - reports/cancer_pc_methods_20261104/README.md
  - reports/cancer_pc_methods_20261104/additional_seed_summary.csv
  - reports/cancer_pc_methods_20261104/development_results.csv
  - reports/cancer_pc_methods_20261104/verification.json
  - reports/cancer_pc_methods_20261104/method_comparison.png
  - reports/cancer_pc_methods_20261104/population_recovery.png
  - reports/cancer_pc_statistics_debug_20261104/README.md
tags:
  - diffusion
  - benchmarks
  - choir
---

# Cancer PC Surrogate Evidence

## Summary

A registered 496-row, 30-feature Gaussian mixture with six unequal populations
reproduces the cancer-inspired singleton-root/closed-edge failure. A matched
one-observation control restores six clusters. Four cases cross the designated
row's noise scale with the separation of two populations, retaining one tree.

## Key Points

- The development hub/disconnected case has graph components 368/64/64,
  a 1/495 root, parent-total covariance fallback, root adjusted edge p-values
  rounding to one, 988 ancestor-blocked edges and one output cluster.
- Changing only the designated row from scale 0.03 to 1 yields a 47/449 root,
  six clusters and ARI 0.998223 on the same random draws and labels.
- Five subsequent seeds all reproduce one versus six clusters. Three retain
  the singleton root. Two have 16/480 roots with significant edges but closed
  sibling gates; p-values 0.876715 and 0.904407.
- For one latter seed the branch-adjusted root statistic changes only from
  316.8763 to 338.0595 in the matched control, while empirical-null inflation
  changes from 14.8439 to 2.35335. Both support records remain below native
  support thresholds, whose enforcement is disabled in this registered method.
- All disconnected runs have three stationary eigenvalues. All root NNLS
  columns are pairwise identical. Native NNLS reports success by cost change;
  full-scale nonconvergence is not reproduced or independently KKT-certified.
- The points above describe code before merge `50455316`. With the
  singleton-peel traversal fix, the development hub/disconnected case keeps its
  1/495 root but passes through it: six clusters, ARI 0.9982. Hub/connected
  moves from four clusters (ARI 0.445980) to five (0.5934); both no-hub cases
  are unchanged. Additional seeds were not rerun.

## Evidence

### Statistical rejection and stopping diagnosis

The follow-up [[cancer-pc-statistics-diagnosis]] traces both raw null-hypothesis
rejection and closed gates on eight native runs. It isolates fitted covariance
and PCA reference mismatch, non-null contamination of stopped-edge calibration,
known row-noise covariance sensitivity, fitted-tree dependence and an exactly
zero-mean split above three symmetric populations. Oracle interventions and
fixed-pair controls diagnose mechanisms without validating a replacement test.

### Default-method comparison and plots

The `cancer_pc` suite now selects these four cases through the standard full
runner. The eleven default methods were evaluated on the development seed and
five additional seeds: 264 rows, comprising 216 successful runs and 48 explicit
continuous-input skips for the two Hamming-based diffusion methods. All 44
development rows match between the standard full runner and the repeated panel.

Across five additional seeds, family-diffusion TBS averages ARI 0 in the
dense-row disconnected case and 0.9922 in its ordinary-row control. Spectral
clustering averages 0.9916/0.9928 and K-means 0.9304/0.9262, with true K=6
supplied to both. Ordinary standardized-Euclidean TBS returns one cluster
throughout. These are registered pipeline comparisons with different distances
and knowledge of K; no parameters were tuned after observing results.

The retained figures show all methods' ARI and cluster-count ranges, plus
development-seed assignments in the two coordinates containing the generating
mean signal. The standard report has 23 tree/UMAP pages. Its common UMAP view
standardizes features and is not the native diffusion-coordinate geometry.
Source/input hashes, every successful ARI/NMI/count and previous family-diffusion
results were independently checked. The suite integration has 24 passing tests.

### Original mechanism reconstruction

The primary report records the model, exploratory parameter choices, frozen
development seed, five subsequent seeds, all four controls, measured scope,
limitations and reproduction commands. The final bundle records exact inputs,
truth/predicted labels, spectra, linkages, gate annotations, calibration replays,
package versions and hashes of the relevant working-tree sources.

The generator enters the existing `method_proof` and `full` suites. Native
input preparation and dispatch are exercised by regression tests. No production
inference or solver implementation was changed for this benchmark.

The simulation begins at the 30-PC representation stage and does not fit
cancer covariance or reproduce counts/preprocessing, 190 populations or rare
component proportions. Its hub has different known row covariance, so the
existing common-sigma conditional chi prototype does not apply. Connected
controls change both graph connectivity and separation. This panel does not
establish selective Type-I error, FDR, exact event certification or a cancer fix.

## Links

- [[cancer-pc-statistics-diagnosis]]
- [[choir-singleton-peel-rerun]]
- [[choir-cell-lines-recovery-plan]]
- [[choir-cell-lines-data]]
- [[conditional-diffusion-monte-carlo-evidence]]
- [[diffusion-topology-branch-rank-isolation]]
