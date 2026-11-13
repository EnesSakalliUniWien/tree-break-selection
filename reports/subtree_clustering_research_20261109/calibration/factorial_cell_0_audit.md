# Independent audit of factorial null cell 0

Date: 2026-11-09, Europe/Vienna. Cell 0 contains 50 independently generated
64-by-12 iid standard-normal datasets. The 16 configurations are paired
on each dataset, giving 800 completed fits. Every root rejection and
every final partition with more than one cluster is a false split for this
single-population truth. Root and partition endpoints are counted
separately. The traces cover replicate 0 in all 16 configurations.

## Canonical dispatcher parity

I regenerated replicate 0 with the factorial seed path
`data_seed(phase=1,panel=3,cell=0,replicate=0)`, giving seed 2328693789.
I used `run_clustering_result` with the first registered parameter grid
entry for each method and the runner's continuous feature space. For
`tbs_family_diffusion` and registered opt-in
`tbs_family_diffusion_markov_ward`, status, complete label arrays and root
sibling p-values matched the factorial's NNLS plus .90 contrast-energy
configurations exactly. Symmetric/average returned three clusters with
root p=2.53770571744617e-17. Markov/Ward returned four clusters with root
p=1.1467092316711193e-35. This establishes the two configured dispatch
paths on the saved null fixture.

## Denominators and paired decisions

Counts below are per 50 inputs for each setting. Finite p counts refer to
recorded root sibling p-values. A missing p-value means the root had no
active numeric sibling p-value, with its exact reason determined from
the trace where available.

| Operator | Linkage | Branch | Rank energy | Root rejects | Finite root p | Final >1 cluster |
| --- | --- | --- | --- | ---: | ---: | ---: |
| Symmetric | Average | Linkage | Off | 35 | 39 | 36 |
| Symmetric | Average | Linkage | On | 38 | 41 | 38 |
| Symmetric | Average | NNLS | Off | 43 | 43 | 43 |
| Symmetric | Average | NNLS | On | 47 | 48 | 48 |
| Symmetric | Ward | Linkage | Off | 4 | 4 | 4 |
| Symmetric | Ward | Linkage | On | 5 | 6 | 5 |
| Symmetric | Ward | NNLS | Off | 40 | 42 | 40 |
| Symmetric | Ward | NNLS | On | 46 | 49 | 47 |
| Markov | Average | Linkage | Off | 27 | 34 | 29 |
| Markov | Average | Linkage | On | 32 | 41 | 34 |
| Markov | Average | NNLS | Off | 32 | 34 | 32 |
| Markov | Average | NNLS | On | 41 | 45 | 44 |
| Markov | Ward | Linkage | Off | 3 | 3 | 3 |
| Markov | Ward | Linkage | On | 3 | 3 | 3 |
| Markov | Ward | NNLS | Off | 41 | 43 | 41 |
| Markov | Ward | NNLS | On | 47 | 50 | 48 |

At nominal root alpha .01, 47/50 has an exact two-sided 95% interval
[.835,.987]. The two rates coincide on this panel, but the data are
paired and cannot be pooled as independent samples. With .90 energy
enabled and operator and linkage held fixed, switching from linkage
lengths to NNLS created 9, 41, 9 and 44 new root rejections in
symmetric/average, symmetric/Ward, Markov/average and Markov/Ward.
None reversed in the opposite direction. Of those, respectively 6,
41, 1 and 44 had no numeric sibling p under linkage lengths.
The corresponding increases in fragmented final partitions were 10,
42, 10 and 45, with no reverse. Some final-partition reversals therefore
occurred without a new root rejection. The trace path can continue
to descendants, so a root rejection is not a complete partition endpoint.

## What changed in the saved traces

On replicate 0, the branch-source switch retains the linkage matrix and
root children within every fixed operator and linkage. At the root,
multiplying raw T by the sibling branch variance factor g reproduced
the same value across linkage and NNLS to the displayed precision for
all four operator/linkage pairs, with df=12 throughout. This empirical
check is consistent with scalar branch rescaling as the dominant change
in raw T on this fixture. The production ridge and rank thresholds
prevent asserting an exact identity in general.

| Operator and linkage | g, linkage → NNLS | Raw T, linkage → NNLS | gT in both, approximately |
| --- | --- | --- | ---: |
| Symmetric average | 2.599 → 1.214 | 56.232 → 120.385 | 146.139 |
| Symmetric Ward | 7.299 → 1.204 | 14.524 → 88.081 | 106.009 |
| Markov average | 2.537 → 1.187 | 45.262 → 96.726 | 114.828 |
| Markov Ward | 7.490 → 1.270 | 34.231 → 201.910 | 256.396 |

- **Root eligibility.** With energy enabled in Markov/Ward, the
  linkage-length child-edge BH p-values were .311 and .0436, both above
  .001. The root trace said `split_prerequisites_closed`, with no sibling
  p and a one-cluster result. Under NNLS, the corresponding values were
  2.20e-15 and 1.01e-17. The root sibling p then became 1.15e-35 and the
  root split. This is a joint change to the fitted branch variance
  normalization, edge eligibility and sibling score. The availability
  of a sibling p-value is part of the observed reversal.
- **Inflation with the root raw score held observationally constant.**
  In symmetric/Ward with NNLS, energy off and on had identical root
  raw T=88.081, df=12 and root edge BH p-values (.000152,.0000239).
  The off setting used family-baseline inflation c=3.99955 and gave
  p=.03727, leaving a one-cluster boundary. The on setting used local
  inflation c=1 and gave p=1.16e-13, producing four clusters. Here the
  root reversal is in the inflation context. The setting can affect
  lower-node projections and evidence even when the root rank is unchanged.
- **Geometry comparisons.** Switching operator and linkage can change
  the split, covariance context, branch fits and eligibility together.
  The paired branch comparison holds operator and linkage fixed. The
  two-method comparison alone cannot separate operator from linkage.

## Scope and source integrity

The manifest reports no changed source during these 800 fits. Compared
with the current checkout, one manifest source differs. I verified the
archived source SHA-256 and independently compared top-level Python ASTs.
Only `runtime_gate` changed. The factorial path does not call it, and the
change converts runtime-summary booleans for JSON serialization.

These large finite-panel rejection rates contradict nominal .01 behavior
for the tested native selected-root and partition procedure on this null
cell. They do not transfer to arbitrary distributions or all settings.
The common-event Haar/Holm proof concerns a different reference and
requires node-wise common Gaussian means and covariance, independent
rows and the actual selection event. Its mixed-family development panel
had no true rejections in 50 inputs at 999 steps. That power finding is a
separate obstacle to promotion. Neither the theorem nor the development
panel rescues the native empirical-inflation procedure or establishes
final-partition control.
