---
title: Diffusion Selected Root Diagnostic
type: analysis
status: reviewed
updated: 2026-11-04
sources:
  - raw/inbox/diffusion-root-cells-certificate-20261102.md
  - raw/inbox/diffusion-root-probability-check-20261102.md
  - raw/inbox/diffusion-positive-self-artifact-20261102.md
  - benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_diagnostic.py
  - benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_kernel.py
  - benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_cells.py
  - benchmarks/diagnostics/calibration/selected/hierarchy/diffusion_root_mass.py
  - tree_break_selection/space_separation/diffusion.py
  - raw/assets/diffusion-root-diagnostic-20261102/final-small-controls/manifest.json
  - raw/assets/diffusion-root-diagnostic-20261102/final-small-controls/roots.csv
  - raw/assets/diffusion-root-diagnostic-20261102/final-small-controls/cells.csv
  - raw/assets/diffusion-root-diagnostic-20261102/final-attribution-first-root/manifest.json
  - raw/assets/diffusion-root-diagnostic-20261102/final-attribution-first-root/roots.csv
  - raw/inbox/diffusion-self-distance-and-region-tightening-20261102.md
  - raw/inbox/diffusion-factor-bound-tightening-20261102.md
  - raw/assets/diffusion-self-distance-20261102/postfix_operator_v2/comparison_summary.csv
  - raw/assets/diffusion-region-tightening-20261102/final-small-controls/roots.csv
  - raw/assets/diffusion-region-tightening-20261102/final-attribution-first-root/roots.csv
  - raw/inbox/conditional-diffusion-statistics-20261104.md
  - raw/assets/conditional-diffusion-20261104/summary.json
tags:
  - diffusion
  - selection
  - roots
  - calibration
---

# Diffusion Selected Root Diagnostic

## Summary

The known-identity diagnostic is implemented with full features, zero branch
lengths, inflation factor 1 and no contrast-energy rank selection. It conditions
on the complete merge history and keeps all unclassified probability mass.
The identity-self-distance artifact is corrected before epsilon selection.
Tighter exact-real bounds certify substantially wider reentry intervals and a
selected interval around the corrected 80-row observation. Conditional p bounds
still remain uninformative. Native floating-point execution is not certified
identical to the exact-real map at every radius.

## Details

### Conditional Monte Carlo root experiment

The subsequent known-identity Gaussian prototype uses the same full-history
radial replay to draw a fixed number of null proposals and form a rank p-value
from accepted radii plus the observation. In 60 nine-row homogeneous Gaussian
datasets it rejected zero roots at alpha 0.01, compared with 20 nominal chi
rejections. Both signal panels also had zero conditional rejections. All 20
strong-signal cases accepted no radii among 999 proposals, forcing p=1.
The experiment exposes poor sampling resolution and leaves root-partition
conditioning and efficient conditional integration open. The rank theorem
assumes correct membership evaluation. Native floating replay has no
exact-real event certificate. [[conditional-diffusion-monte-carlo-evidence]]
preserves the derivation, Wolfram checks, paired outputs and finite-panel limits.

### Conditional reference law

For independent Gaussian rows with known identity covariance and root contrast
eta satisfying eta.T E[X]=0, set v=X.T eta and H=eta.T eta. The exact path is

\[
X(r)=X-\eta v^T/H+\frac{\eta v^T}{\sqrt{H\,v^Tv}}r,
\qquad r_{obs}=\sqrt{v^Tv/H}.
\]

Conditional on the nuisance offset, contrast direction and full merge history
M, the radius has chi_p density restricted to S_M. Here p is the input feature
count; native diffusion-coordinate truncation remains part of construction.
The diagnostic separately encloses the observed radius and every operator
quantity. It never substitutes the rounded replay path for this exact path.

For certified selected mass A above r_obs and B below, and unresolved mass U
above and V below,

\[
\frac{A}{A+B+V}\le p_M\le\frac{A+U}{A+U+B}.
\]

Arithmetic and endpoint uncertainty widen these bounds. Without a positive
selected-mass lower bound and certified observed-event membership, the gate
uses [0,1]. It rejects only if the upper bound is at most alpha; otherwise a
straddling bound retains one cluster without passthrough. This is root-only
validation, not a descendant traversal or multiplicity guarantee.

### What interval certification currently proves

The kernel encloser supports the native dense route, presently 9 through 2,000
rows, subject to proving all required adaptive branches. The upper size bound
is an applicability limit, not a performance claim. Spectral certificates
require an isolated leading mode, a separated retained cutoff, and positive
retained eigenvalues. Unproved branches and arithmetic failures are unresolved.

For seed 20261105, n=9, p=1, windows around radius ratios 1, 1.1 and 1.2 are
respectively selected, nonselected and selected. On identical operator boxes,
the largest tested certified selected half-width increases from 1e-12 with the
old Gram bound, through 1e-7 with a global diffusion-factor bound, to 1e-5
with the retained-spectrum factor bound. These are tested widths, not recovered
selection boundaries. All gaps and the infinite tail retain unknown mass.

The improvements come from exact pair-distance quadratics (shared within-child
motion cancels), proving only the required neighbor support sets, and bounding
distances as ||H(e_i-e_j)|| for H=S^t P. Projecting the power telescope onto
retained spectra avoids including the removed leading mode in every norm.
Validated expansion accepts a larger interval only after a whole-interval proof;
bounded subdivision refines unresolved finite gaps.

Final scalar panels use half-width 1e-5, four expansion attempts and subdivision
depth 4. For the first control, selected-mass lower bound grows from about
8.81e-14 to 1.05858e-6; p bounds are approximately [7.16085e-7,0.999972].
The third small control still lacks certified observed-event membership.
The corrected first 80-by-24 attribution root has a certified selected interval
and selected-mass lower bound about 4.88479e-14. All four gates remain
unresolved and retain one cluster. This is not successful calibration.

### Native self-distance defect and its correction

Before correction, the first original 80-by-24 selected matrix had 29 slightly
positive self distances. Native epsilon took the median of all positive scaled-distance
squares before the weight diagonal is zeroed. Consequently these artifacts
enter epsilon. The reconstructed radial observation has 28 such entries.

| Native input | Positive self entries | Native epsilon | Off-diagonal median |
| --- | ---: | ---: | ---: |
| Original selected matrix | 29 | 0.8588652016711084 | 0.8639359788816151 |
| Reconstructed radial observation | 28 | 0.859336163375187 | 0.8639359788816152 |

In the pre-correction investigation, a floating eigensolver/linkage proposal
for the zero-self operator first changed the history at merge index 5.
The corresponding interval calculation strengthened this:
a positive-width cell containing the exact observed radius is certified
nonselected for the native history. At step 5, competitor pair [22,66] beats
the chosen [9,72], with an upper distance-difference bound of about -5.55e-7.
This excludes the native history without identifying the entire alternative
history. That earlier observed event was therefore unsuitable for the exact-real
law. Self-distance counts depend on numerical execution; provenance accompanies
the fixture.

The production constructor now zeros its identity diagonal after kernel fitting
and before epsilon selection/kernel computation. It preserves every off-diagonal
entry, including stored zero distances with unit kernel affinity. The corrected
original/reconstructed first-root epsilon values are respectively
0.8639359788816151 and 0.8639359788816152. Both have zero positive self entries.
The corrected observed history now has a positive selected-cell certificate.

Genuine before/after observations cover 50 Gaussian replicates, three small
controls and two continuous signal cases. Eight of 50 Gaussian complete
histories change ordering, while no rooted clade sets or root partitions change.
The other five histories remain unchanged. All original/reconstructed pairs
match within each version. Known-identity root statistics are exactly unchanged,
and descriptive unselected chi-square tails still reject all 50 Gaussian roots
at .01. The self-distance correction therefore does not remove selection bias.

The exact-real certificate explicitly uses zero self distances and rational
bandwidth exponents. It is not a proof of bitwise NumPy/NumExpr execution.
Only self-distance handling changed in production; thresholds and traversal
were not changed. KDE fitting still uses native floating arithmetic.

### Reproduction

```bash
uv run --no-sync python -m benchmarks.diagnostics.calibration.selected.hierarchy.diffusion_root_diagnostic --replicates 3 --observed-half-width 1e-5 --expansion-steps 4 --subdivision-depth 4 --output-dir /tmp/diffusion-root-diagnostic
```

The runner writes separate root decisions, interval cells and point replays,
with package versions, input/source hashes and unresolved-case counts. Point
observations explicitly carry `interval_certificate=False`.

## Evidence

[[diffusion-root-diagnostic-evidence]] indexes the derivations, Wolfram capture,
five implementation modules and five test modules. The native artifact note
includes the 80-row source-level diagnosis and reproduction commands.

Final panel artifacts are under
`raw/assets/diffusion-root-diagnostic-20261102/final-small-controls/` and
`final-attribution-first-root/`. Earlier sibling directories preserve intermediate
implementation runs with their own source hashes; they are not the final evidence.
Those directories record the earlier implementation stage. Current corrected
topology evidence is `raw/assets/diffusion-self-distance-20261102/postfix_operator_v2/`;
current tighter panels are the `final-*` directories under
`raw/assets/diffusion-region-tightening-20261102/`. The subsequent complete
ordered suite passes 977 tests. [[diffusion-self-distance-repair-evidence]] links
the correction, isolated comparisons, new derivations and limitations.
Ruff lint, changed-file formatting, audits, generated-file checks and 228-page
wiki lint pass. Repository-wide formatting still flags the unchanged earlier
raw reentry probe; it was preserved as evidence.

## Links

- [[diffusion-root-diagnostic-evidence]]
- [[diffusion-self-distance-repair-evidence]]
- [[diffusion-selected-law-open-question-answers]]
- [[diffusion-selection-covariance-attribution]]
- [[analytic-selected-root-laws]]

## Open Questions

- Reconcile remaining numerical uncertainty with the declared selection model
  before claiming bitwise native-map validity at all radial points.
- Improve interval bounds enough to cover substantial conditional mass and
  handle branch boundaries, ties, tangencies and tails efficiently.
- Validate calibrated decisions under the resulting declared constructor;
  unresolved outcomes cannot be counted as successful nonrejections.
- Derive estimated-covariance, fitted projection/inflation and descendant
  traversal laws only after the known-identity root construction is resolved.
