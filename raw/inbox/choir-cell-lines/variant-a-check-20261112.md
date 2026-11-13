# Independent check of variant A of the draft definition of k

Check on 2026-11-12, Europe/Vienna. Local folder: `results/analyses/variant_a_check_20261112` (git-ignored); the
brief, the checker's result summaries, its two central scripts and the own verification scripts with their
outputs are in `variant-a-check-scripts-20261112.py.txt`. Follow-up of `k-split-direction-level-20261112.md`.
No production code was changed.

## Setup

The brief with six questions was written before the check. The check was done separately from the earlier work
on variant A; each finding used below was then recomputed or rerun with separate code, except where marked
"checker's measurement". The checker did not produce a report file; this capture is the record.

## Findings

### Defect: an error that stops the fit (confirmed by an own call)

When a feature is constant in the basis half and the basis half has more than p + 2 cells, the pooled within
covariance is singular and the draft raises `LinAlgError`; `apply_best_flat_cut_reference` does not catch it, so
the whole fit stops. The current rule returns p = 1.0 on the same input. It did not occur in any run so far (the
benchmark data and the 30 PCs have no constant feature inside a node). Variant A has to fail closed here.

### Limitation: the result depends on the order of the rows (confirmed by an own run)

The halves alternate in the order of the rows, so the same cells in another order give other halves.

CHOIR cells, Markov/Ward tree, variant A:

| Order of the rows | Found clusters | ARI | Lines recovered | Accepted by BH |
| --- | ---: | ---: | ---: | ---: |
| Data order (the run of the same day) | 218 | 0.8311 | 120 | 151 |
| Random order 1 | 225 | 0.8397 | 125 | 156 |
| Random order 2 | 216 | 0.8321 | 123 | 149 |
| Random order 3 | 205 | 0.8236 | 118 | 144 |
| Random order of the own run | 196 | 0.8043 | 114 | 139 |
| Current rule (no dependence on the order) | 230 | 0.8219 | 116 | 134 |

- Over the first four orders, 131 nodes are accepted in every order and 39 in only some; of 191 nodes with
  p <= 0.01 in some order, 25 have p = 1 in another; the partitions of two orders have an ARI of 0.9585 to 0.9847.
- The earlier statement "ARI 0.8311 against 0.8219" was one draw: variant A is above the current rule in 4 of 5
  orders and below it in 1.
- Benchmark, 14 cases, 5 random orders: found clusters and ARI are unchanged in 28 of 28 rows in every order;
  the set of accepted splits changes in 0, 4, 4, 2 and 3 of 28 rows.

### Limitation: one population that is not Gaussian (checker's measurement, summary read, not rerun)

256 and 496 cells in 30 features, 20 seeds, 4 operator and linkage cells, 160 fits per shape and rule:

| Shape of the one population | Fits with more than 1 cluster: variant A | Current rule |
| --- | ---: | ---: |
| Skewed features | 3 of 160 | 1 of 160 |
| Heavy tails | 0 of 160 | 0 of 160 |
| Uniform | 1 of 160 | 2 of 160 |
| Closed curve of radius 3 in two features, noise in the other 28 | 143 of 160 | 158 of 160 |

At 496 skewed cells variant A has p <= 0.01 at 24 of 2,200 tested nodes (1.09%; 0.35% with one Gaussian
population). The closed curve is split by both rules; this is a property of the Gaussian reference, not of k.

### No problem found

- Code against the definition: the four points are implemented as written; sides of 2 or 3 cells give the same
  direction as the pooled formula; an infinite `rho` gives p = 0.0005; test cells are rounded down to the table's
  grid.
- CHOIR cells, composition of the accepted nodes (recomputed with own code): variant A accepts 151 splits, 3 at
  nodes with at least 90% of the cells from one line; the current rule 134, of them 8. Of the 27 splits only
  variant A accepts, 1 is inside one line; of the 10 only the current rule accepts, 6.
- One Gaussian population, recorded level data: p <= 0.1, 0.05 and 0.01 at 1.87%, 1.04% and 0.35% of 31,482
  tested nodes; every class of test cells is at or below the nominal share.
- The current rule's run reproduces exactly, and Benjamini-Hochberg sees the same 1,838 nodes under both rules.

### Limitation of the production table at k = 1 (confirmed by an own simulation)

The exact best cut of one direction, 20,000 Gaussian datasets per n, share with reference p <= 0.01: n = 15:
1.65%; n = 30: 1.03%; n = 75: 1.12%; n = 100: 1.89%; n = 150: 0.66%; n = 400: 1.40%. The checker's 100,000
datasets per n give 1.62%, 1.84% and 1.41% at n = 15, 100 and 400 and 0.66% to 1.20% elsewhere. The table has
1,999 draws per cell, so its 99% points carry this error. It concerns the opt-in reference in general; variant
A's splits lie far below the best cut, so its measured level is hardly affected.

### Correction of one of the checker's findings

The checker reported the benchmark mean "0.8302 against 0.6950" as mislabelled. Recomputed: it is the mean over
the 9 continuous cases of 11 November and leaves out the null case added on 12 November
(`gauss_null_256n_30f_continuous`), as the capture of that run says. Over all 10 continuous cases the means are
0.8472 and 0.7255; over all 14 cases 0.8114 and 0.7245. No number was wrong; messages that said "continuous
suite" without "9 cases" were imprecise.

## Reading

- Variant A's code and its level under one Gaussian population stand.
- Two things have to change before it is used: it must fail closed where the direction cannot be computed, and
  the halves must not depend on the order of the rows, or the result must be reported with that spread.
- On the CHOIR cells variant A and the current rule are not distinguishable: 196 to 225 clusters and ARI 0.8043
  to 0.8397 against 230 and 0.8219. What remains in favour of variant A is the benchmark (three cases the
  current rule loses) and that it accepts fewer splits inside one cell line (3 against 8).

## Not checked

The non-Gaussian populations were not rerun with own code. More than 5 orders on the CHOIR cells; the
one-population grid under other orders; other trees and datasets; data with features that are constant inside a
node; a way of forming the halves that does not depend on the order.
