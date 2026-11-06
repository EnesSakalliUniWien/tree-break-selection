# Cancer-inspired PC benchmark: root failure and matched controls

Date: 2026-11-04. Implemented benchmark generator and four registered method-proof
cases; no production clustering, inference, solver, or threshold changes.

## Result

A 496-row, 30-feature synthetic case reproduces three disconnected diffusion
components, a densely connected singleton at the root of one average-linkage
tree, parent-total covariance fallback, two closed root edge gates, and a
one-cluster result that prevents exploration of six generating populations.
Restoring ordinary noise to just that observation produces six clusters.

The five subsequent seeds all retain the one-cluster versus six-cluster
difference. The exact singleton-root mechanism occurs in three of those five;
the other two stop at the root sibling gate after significant root edge tests.
This is a mechanism benchmark with known generating structure, not a fitted
replica of the cancer experiment or an inference-calibration result.

## Frozen model

For independent standard normal 30-vectors Z_i, generate

    X_i = mu[z_i] + s_i Z_i.

Population sizes are [256, 64, 32, 16, 64, 64]. The first four means have first
coordinate [0, 6, 12, 18]. The last two have second coordinate [-a, a]; all
other coordinates of every mean are zero. Ordinary rows have s_i=1. The
designated first observation S0000 has s_0=0.03 in the hub case, or 1 in the
matched control. Both controls use the same random draws and population labels.
Use a=20 for the disconnected case and a=3 for the connected control. Every
population still has a different mean in both configurations.

The generator starts at a continuous representation corresponding to the
30-PC stage. It does not generate gene counts, perform normalization/HVG
selection/PCA, or estimate cancer-data covariance. The 64-row isolated groups
are deliberately enlarged relative to the original 250/185-row components
among 48,879 cells so that native k=15 neighbor construction remains meaningful.
There are six populations, versus 190 annotated cancer cell lines. With
s_0=0.03 the known row covariances differ: Cov(X_0)=0.0009 I, Cov(X_i)=I for i>0.

Native `tbs_family_diffusion` uses its registered parameters, k=15, diffusion
time 3, symmetric coordinates, average linkage, and NNLS targets in original
feature coordinates. Edge alpha is 0.001 and sibling alpha is 0.01. One tree
contains all observations. The original continuous data are never replaced by
diffusion coordinates for the tests. The runner's current working-tree
empirical-null fallback is included in the captured source provenance.

## Development and held-out-seed results

One development seed, 20261104, was used to choose a reproducing model.
An initial grid with main sizes [72,28,16,12], isolated sizes 24 or 48,
main separation 2 or 4, and hub scale 0.03 or 1 did not produce the singleton
root. Increasing the main sizes to [256,64,32,16], with isolated sizes 64 or
128 and separation 4 or 6, did. The smaller 496-row model with separation 6
was frozen before examining seeds 20261105 through 20261109. There was no
search for favorable p-values or modification of inference thresholds.

| Development case | Graph components | Root child sizes | TBS clusters | ARI |
| --- | --- | --- | ---: | ---: |
| Hub, disconnected | 368 / 64 / 64 | 1 / 495 | 1 | 0 |
| Ordinary row, disconnected | 368 / 64 / 64 | 47 / 449 | 6 | 0.998223 |
| Hub, connected | 496 | 16 / 480 | 4 | 0.445980 |
| Ordinary row, connected | 496 | 16 / 480 | 6 | 0.697428 |

| Five additional seeds | Singleton-root runs | TBS cluster counts | ARI range |
| --- | ---: | --- | --- |
| Hub, disconnected | 3/5 | 1, 1, 1, 1, 1 | 0 |
| Ordinary row, disconnected | 0/5 | 6, 6, 6, 6, 6 | 0.986397–1 |
| Hub, connected | 0/5 | 4, 4, 4, 4, 4 | 0.428987–0.447305 |
| Ordinary row, connected | 0/5 | 6, 4, 4, 5, 6 | 0.428987–0.661441 |

Connectivity and separation both change when a changes. The connected controls
therefore do not isolate a graph-only effect or establish that reconnecting a
real graph is a remedy. The one-row control also changes that observation's
contribution to statistics; the benchmark isolates a data intervention, not
the individual causal contribution of each downstream operation.

## What the benchmark explains

### Singleton root and blocked descendants

In the development hub/disconnected case S0000 has 256 undirected graph
neighbors, 255 from its own population. It is not an isolated graph vertex.
Yet symmetric diffusion plus average linkage leaves it alone at the root.
The mean standardized squared distance across that root is 1.0102 versus
2.0040 within the larger child. Thus the large diffusion-coordinate separation
does not imply a large original-feature mean contrast.

The root has a child smaller than the covariance guard of eight leaves. Its
resolved covariance equals parent-total covariance; mean diagonal 5.1591
contains population mean structure. Both root adjusted edge p-values round
to 1, and only two edges are tested. There are 988 ancestor-blocked edges and
one live visited node. Returning one cluster here does not establish that
the dataset is homogeneous.

The original cancer singleton had 88 neighbors, all from ASPC1, root adjusted
p-values 0.9167/0.7246, and parent-total diagonal covariance. The small surrogate
uses a full 30-by-30 covariance and has different cross/within distances and
population proportions. These are explicit limits to the match.

### Sibling calibration supplies a second failure route

Seeds 20261106 and 20261109 have a 16/480 root in both disconnected controls.
In the hub cases both root edge tests reject, but sibling p-values are
0.876715 and 0.904407. Traversal still visits only the root.

For seed 20261106, the hub root's branch-adjusted statistic is 316.8763 and
the fitted inflation factor is 14.8439, producing a reference statistic
21.3473 on 30 degrees of freedom and p=0.876715. Its matched control has
statistic 338.0595 and inflation 2.35335, producing 143.6501 and p=8.844e-17.
The hub uses the family-baseline fallback; the control uses the local kernel.
The empirical-null support metadata records below-threshold support in both
cases, while native enforcement is disabled. Recomputed calibration decisions
match the actual gate p-values. This is descriptive behavior of the present
implementation, not evidence that either selected reference law is valid.

### Repeated eigenvalues and NNLS identifiability

All 12 disconnected runs have three eigenvalues numerically equal to one.
The connected controls have one. This recreates the repeated stationary
subspace that violates the current exact-real certifier's isolated-leading-mode
assumption; this experiment does not itself execute a new event certificate.

The two root-edge NNLS design columns are exactly identical in every run.
Only their sum is identifiable from leaf-pair distances. This structural fact
also holds in the cancer trees and does not alone imply solver failure.

All 24 native fits report success using their cost-change stopping criterion.
The development hub fit uses 50,000 of 122,760 pairs and samples 208 crossings
of its singleton root. Its reported optimality is 2.5345, so this report does
not claim an independently verified KKT certificate. A separate sampling-only
control scales the budget to 507 pairs to match the cancer pairs-per-cell
ratio; at seed zero it samples no root crossing. No TBS result uses that
reduced-budget fit. Pair-sampling fragility is exposed, but full-scale NNLS
iteration-limit failures and dense-memory pressure are not reproduced here.

## Conditional selective inference scope

The hub cases violate the existing prototype's common known spherical
covariance assumption. Their heteroskedastic covariance is known by
construction, making them useful for a future correctly derived path, but
the current common-sigma chi law must not be applied to them. The nohub cases
have identity noise covariance; their selected root contrasts have nonzero
population means. This signal panel estimates neither Type-I error nor FDR.
Full-history Monte Carlo acceptance, descendant selective testing and exact
event certification were not rerun in this experiment.

## Reproduction and evidence

The generator is `benchmarks/shared/generators/cancer_pc_surrogate.py`.
The four `cancer_pc_*` recipes are registered in the existing `method_proof`
and `full` suites. The curated nine-case `continuous` suite is unchanged.

Load the principal case through the usual entrypoint:

```python
from benchmarks.shared.cases import get_test_cases_by_suite
from benchmarks.shared.generators.generate_case_data import generate_case_data

case = next(c for c in get_test_cases_by_suite("method_proof")
            if c["name"] == "cancer_pc_hub_disconnected")
data, truth, original, metadata = generate_case_data(case)
```

Replay all 24 runs into a new directory from the repository root:

```bash
PYTHONPATH=. uv run --no-sync --offline python \
  raw/inbox/cancer-pc-surrogate-experiment-20261104.py.txt \
  /tmp/cancer-pc-surrogate-replay
```

`raw/assets/cancer-pc-surrogate-20261104-final/` contains the authoritative
manifest with package versions/source hashes, 24 result records, a compact
summary CSV, and four compressed development inputs with truth, predicted
labels, linkage, spectrum and components. Full development annotations and
root annotations/calibration decisions for every seed preserve both stopping
routes. `raw/assets/cancer-pc-surrogate-20261104/` is the earlier panel before
adding explicit sibling-calibration replay and formatting the implementation.

`tests/core/test_cancer_pc_surrogate.py` checks registration, native input
preparation, deterministic paired controls, finite input validation, three
components/repeated leading modes, and the native one-cluster/six-cluster
regression. The tests first failed because the generator and cases were absent.
The existing 18 case-contract tests passed before edits; the combined focused
suite subsequently passed 23 tests. Final lint and artifact verification are
recorded with the evidence bundle.
