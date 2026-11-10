# Root and immediate-child calibration follow-up

Date: 2026-11-09, Europe/Vienna. Scope: saved replay, independent known-noise
Gaussian panels, a paired chain budget, and isolated A/B/C cut diagnostics.
The runnable protocols, scripts, full chains, per-cut records, summaries and
manifests are under `results/analyses/root_validation_followup_20261109/calibration/`.
That generated directory is git-ignored. The companion source archive retains
the executable code and locked protocols.

## Saved replay and fresh calibration

All 24 saved 999-step root chains reproduced their stored radii, random anchor,
rank, moves, status and decision exactly. The first replicate of all seven
saved root-plus-child panels reproduced the stored per-node summaries
(observed radii checked at relative tolerance 1e-13). All 45 tests from the
three existing sampler test files passed. Verification also checked every
complete fresh rank and anchor and all six recorded source hashes.

All 1,664 fresh datasets ran at 119 steps, with the observed statistic at an
independent random chain position. There were 4,441 moved chains, five stalled
chains and one unsupported observed tree with a numerical tie. Stalled and
unsupported attempts remain in dataset denominators. No seed or budget was
replaced. The root-only cutoff uses the same root chain as the family. The
family tests the root and its internal children, each at 0.05/3. False
rejections use the true contrast means for the selected groups.

| Panel | n, p | Root at .05 | Root at .01 | Any false family rejection | Any true family rejection | Exact top split |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Known-noise small null | 9, 3 | 38/1200 | 1/1200 | 19/1200 | | |
| Known-noise wide null | 12, 8 | 6/120 | 1/120 | 3/120 | | |
| Known-noise larger null | 24, 5 | 6/120 | 2/120 | 7/120 | | |
| Balanced shift 4 | 24, 5 | 14/48 | 3/48 | 1/48 | 7/48 | 10/48 |
| Imbalanced shift 4 | 24, 5 | 12/48 | 2/48 | 0/48 | 7/48 | 7/48 |
| Strong four-group hierarchy | 32, 5 | 31/32 | 6/32 | 0/32 | 16/32 | 32/32 |
| Fitted-scale null sensitivity | 9, 3 | 2/96 | 0/96 | 1/96 | | |

In the primary 1,200-dataset panel, exact 95% two-sided intervals were
[.02250,.04321] for the root at .05, [.0000211,.004634] for the root at .01,
and [.00956,.02462] for the family. One-sided upper limits were .04130,
.00395 and .02315, respectively. Both .05 endpoints met the fixed .065
practical compatibility threshold and their intervals were below nominal
in this panel. The wider and larger secondary family upper limits were
.06334 and .10677, respectively. These estimates concern the tested model,
geometry and budget. The fitted-scale panel supplies a sensitivity result.

In the four-group hierarchy, the exact top split recovered the intended two
halves on all 32 inputs. Immediate-child signal rejections were 2/28 and 2/32
by signature position, and both sides of the selected cut were pure on only
7/28 and 2/32 of those true-signal positions. Median root moves ranged from 12
to 32 across panels. Median total family runtime ranged from 0.85 to 1.45
seconds on this machine with four worker processes and single-thread BLAS.

## Paired chain-budget check

For the first 20 balanced signal inputs, moving from 119 to 999 steps raised
root rejection at .05 from 6/20 to 12/20 (six gains, no losses) and at .01
from 1/20 to 8/20 (seven gains, no losses). Any true family rejection rose
from 2/20 to 8/20 (seven gains, one loss). The same single family had a false
rejection at both budgets. Median root moves rose from 12 to 84, and median
family runtime from 1.42 to 11.65 seconds.

Among the first 30 primary null inputs, root rejection at .05 remained 1/30
on the same input. Family false rejection was 1/30 at each budget on different
inputs. There were no root rejections at .01, no unresolved pairs, and the
median root moves rose from 29.5 to 263.5. Median family runtime rose from
1.09 to 9.39 seconds. All paired input hashes, seeds and selected splits
matched. These 50 pairs quantify a budget effect on these inputs.

## A/B/C cut diagnostics

There were 1,000 datasets per panel and two cuts per dataset. The fixed labels
were assigned before data generation. The other cut maximized the Fisher ratio
over thresholds along a single stated input coordinate. This search is a
subset of all flat cuts. A denotes the isolated projected-Wald kernel and
does not include production branch times, empirical-null inflation or
traversal. B/C use all supplied coordinates and their ceiling-3 references.
The Gaussian-shape C comparison is retained separately.

At .01, on fixed unit-Gaussian labels at n=60,k=6, the known-spherical A
kernel rejected 9/1000, pooled-covariance chi-square A rejected 25/1000,
and its 90%-energy version rejected 31/1000. Exact Hotelling F rejected
5/1000. Under a paired condition-number-100 invertible transform the
misspecified spherical A rejected 65/1000, while the fixed-label B/C/F
values matched their paired originals. The heavy-tail spherical A rejected
23/1000. On the searched Gaussian thresholds the unadjusted F and A
variants rejected 1000/1000, illustrating selection directly.

At n=32,k=20, fixed unit-Gaussian labels gave 12/1000 for known-spherical
A, 763/1000 for estimated pooled chi-square A, and 8/1000 for F. The searched
threshold gave 818/1000 for B and 955/1000 for the unadjusted F reference.
C clips to 1 at every attainable statistic in this setting. The bound R<1
gives u<4.793221 for ceiling 3, where the monotone envelope still exceeds
4.02 million before clipping. The Gaussian-shape C envelope also exceeds
one at its corresponding upper endpoint.

Refitted count normalization, log1p, high-variance feature selection,
scaling and PCA were run on simulated counts with a shared underlying
composition. B/C rejected 0/1000 for both cuts at constant depth. When two
depths were assigned in random order, median absolute PC1/log-depth
correlation was .929 and the fixed-label B/C rejection counts remained
zero. The searched cut gave 977/1000 for B and 143/1000 for C. When the
fixed labels coincided with the two depth levels, B/C gave 771/1000 and
60/1000 for fixed labels and 966/1000 and 163/1000 for the searched cut.
These numbers measure sensitivity under a shared-composition biological null.
Unequal depth can alter the transformed feature distribution.

The chain panels provide evidence for root and immediate-child tests on raw
independent Gaussian features with known spherical noise. Count preprocessing,
estimated covariance, arbitrary descendants, pass-through, full traversal and
final-partition error control remain outside that validation.

For fixed groups under known spherical Gaussian noise, the null concerns their
true mean contrast. With a unit row contrast e, the noise in e'X is independent
of (I-ee')X. Under e'E[X]=0, its radius in sigma units is chi-distributed and
independent of its direction. This remains true when row means differ within
the residual space. Conditioning on the selected ancestor path restricts that
radial law. This zero-contrast null differs from the common-mean null of the
separately studied unknown-covariance Haar route.
