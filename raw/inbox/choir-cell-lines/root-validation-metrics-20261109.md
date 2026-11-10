# Matched partition agreement and robustness

## CHOIR cell lines

All eight saved fits cover the same 48,879 cell IDs in the same order and all 190 cell-line labels. Cluster counts, legacy 90/90 counts, ARI and NMI against cell lines, and ARI against published CHOIR_1 reproduce. Both archived rule-out assignments also match their previous archive cell by cell. The metrics use a one-to-one maximum-total-F1 match to truth, with zero precision, recall and F1 for unmatched lines. Macro metrics weight each line equally. The rare-line set is the 48 smallest lines (56–145 cells), with lexical ties. Legacy 90/90 uses each line's largest-overlap cluster and that cluster's maximum cell-line purity.

| Tree | Test | Clusters (signed error) | ARI | AMI | Homogeneity | Completeness | Macro F1 | Rare F1 | 90/90 (rare) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Symmetric/average | A | 75 (−115) | .028 | .533 | .385 | .949 | .184 | .089 | 27 (3) |
| Symmetric/average | B | 64 (−126) | .028 | .535 | .386 | .956 | .194 | .128 | 29 (4) |
| Symmetric/average | C | 57 (−133) | .028 | .536 | .384 | .962 | .189 | .104 | 28 (3) |
| Symmetric/average | rule-out | 69 (−121) | .018 | .480 | .337 | .943 | .194 | .130 | 29 (4) |
| Markov/Ward | A | 954 (+764) | .549 | .831 | .956 | .786 | .617 | .701 | 41 (16) |
| Markov/Ward | B | 217 (+27) | .833 | .935 | .931 | .948 | .800 | .715 | 119 (25) |
| Markov/Ward | C | 173 (−17) | .815 | .933 | .919 | .957 | .746 | .557 | 111 (20) |
| Markov/Ward | rule-out | 664 (+474) | .765 | .884 | .954 | .853 | .741 | .804 | 80 (23) |

B on Markov/Ward reduces the share of same-line pairs assigned apart from .572 under A to .085. The share of same-cluster pairs from different lines is .233 for B, .227 for A, .270 for C, and .168 for rule-out. C's smaller absolute count error (17 vs 27 for B) coincides with worse macro and rare F1. For the smallest lines B has higher 90/90 recovery (25 vs 23) and lower matched mean precision and recall (.724/.735 vs .890/.771) than rule-out. Symmetric/average merges different cell lines heavily in all four arms (.978–.983 of within-cluster pairs). The split and merge metrics describe label agreement.

## Relationships among saved splits

B and C are evaluated on the same 90 symmetric/average and 1,838 Markov/Ward nodes, with identical k and rho. At raw .01, Markov/Ward has 105 both, 16 B only, 1 C only and 1,716 neither. Symmetric/average has 10 both, 4 B only and 76 neither. On Markov/Ward, no n≤31 node passes, 14 B and 1 C pass among 867 with n=32–159, and 42 B and 43 C pass among 103 with n≥1,000. Smaller-child share, depth and k strata are in `choir_split_strata.csv`. Node size, k, balance and depth co-vary.

The 167 matched Markov/Ward spectrum/distance rows include 13 with different dominant lines and ≥.9 purity on both sides, 20 with a parent ≥.9 pure in one line, and 134 mixed. Raw passes are 11/13 for both, 0/20 for both, and 64/134 B or 65/134 C for mixed nodes. The median test-space centroid-distance tie share is .0051 for the 75 both-pass nodes with distance records versus .0586 for the 91 neither-pass nodes. For invariant-space ties the corresponding medians are .0294 versus .0781. The 10 nodes with top spectral multiplicity 2 yield two passes in each arm, compared with 73 B and 74 C out of 157 with multiplicity 1. This subset was selected by the earlier study's node size requirements.

The previous A annotation partition agrees with the audit's A partition (cell-aligned ARI 1). Its tested-node set matches B/C's eligibility set. On Markov/Ward, the median archived A sibling null weight is the positive floating-point floor (4.94e−324) among B+C passes and about 10^−22 among neither-pass nodes. The respective median A inflation factors are 331.6 and 9.49. The tested-node weights are selected through the edge gate and are heavily floor-censored. They describe the relationship among archived quantities for eligible tests. The earlier reference-weight study treats the reference population separately. These are nested nodes, and the joined labels and scores support descriptive comparisons.

## Fresh fitted partition panel

The frozen panel pairs all four arms on the same 4 seeds, 2 trees, 7 regimes and original/jittered inputs. It rebuilds feature-family diffusion, linkage, branch-time fitting, and traversal. The table uses original inputs. Entries are means over four seeds and follow arm order A / B / C / production rule-out. The full tables include precision, recall, signed count error, split/merge pair shares, runtime, and paired differences from A.

| Regime | Tree | Mean clusters | Mean ARI | Mean macro F1 |
|---|---|---|---|---|
| base | Symmetric/average | 4.00 / 4.00 / 4.00 / 4.00 | 0.964 / 0.964 / 0.964 / 0.964 | 0.986 / 0.986 / 0.986 / 0.986 |
| base | Markov/Ward | 8.00 / 4.00 / 4.00 / 4.00 | 0.720 / 0.961 / 0.961 / 0.961 | 0.793 / 0.985 / 0.985 / 0.985 |
| larger_n | Symmetric/average | 5.50 / 4.00 / 4.00 / 4.75 | 0.913 / 0.961 / 0.961 / 0.865 | 0.952 / 0.985 / 0.985 / 0.885 |
| larger_n | Markov/Ward | 18.75 / 4.00 / 4.00 / 3.00 | 0.344 / 0.953 / 0.953 / 0.651 | 0.522 / 0.982 / 0.982 / 0.682 |
| wider_p | Symmetric/average | 4.00 / 3.50 / 3.25 / 2.25 | 0.872 / 0.793 / 0.731 / 0.394 | 0.949 / 0.810 / 0.732 / 0.454 |
| wider_p | Markov/Ward | 2.75 / 3.50 / 3.25 / 2.00 | 0.609 / 0.783 / 0.739 / 0.329 | 0.567 / 0.809 / 0.731 / 0.375 |
| imbalanced | Symmetric/average | 5.25 / 3.00 / 2.75 / 3.50 | 0.808 / 0.709 / 0.703 / 0.752 | 0.922 / 0.698 / 0.635 / 0.817 |
| imbalanced | Markov/Ward | 11.50 / 2.25 / 2.25 / 3.25 | 0.332 / 0.616 / 0.616 / 0.722 | 0.717 / 0.510 / 0.510 / 0.780 |
| lower_gap | Symmetric/average | 9.75 / 1.00 / 1.00 / 1.75 | 0.276 / 0.000 / 0.000 / 0.083 | 0.560 / 0.100 / 0.100 / 0.197 |
| lower_gap | Markov/Ward | 11.50 / 1.00 / 1.00 / 1.00 | 0.250 / 0.000 / 0.000 / 0.000 | 0.527 / 0.100 / 0.100 / 0.100 |
| anisotropic | Symmetric/average | 8.25 / 1.75 / 1.75 / 3.00 | 0.642 / 0.224 / 0.224 / 0.489 | 0.783 / 0.315 / 0.315 / 0.581 |
| anisotropic | Markov/Ward | 13.00 / 1.75 / 1.75 / 2.75 | 0.391 / 0.226 / 0.226 / 0.534 | 0.583 / 0.316 / 0.316 / 0.580 |
| one_population | Symmetric/average | 8.25 / 1.00 / 1.00 / 2.25 | 0.250 / 1.000 / 1.000 / 0.750 | 0.478 / 1.000 / 1.000 / 0.870 |
| one_population | Markov/Ward | 12.75 / 1.00 / 1.00 / 4.75 | 0.000 / 1.000 / 1.000 / 0.250 | 0.209 / 1.000 / 1.000 / 0.659 |

- Increasing n from 240 to 480 in the balanced gap-5 regime changes the Markov/Ward A mean from 8.0 to 18.75 clusters, while B and C recover four on all seeds and mean ARI .953 at n=480. At p=72, symmetric/average A recovers four clusters on average with ARI .872, while B and C under-count with ARI .793 and .731. These contrasts hold only for this design.
- At 55/25/15/5 imbalance, Markov/Ward B/C raise ARI from .332 to .616 and reduce macro F1 from .717 to .510. The 12-cell component is unmatched on all four original seeds under B/C (rare F1 0). Its mean F1 is .935 under A and .730 under rule-out. At gap 2.5, B/C return one cluster on both trees and every original and perturbed input. With anisotropic covariance they average 1.75 clusters on each tree and original-input ARI .224–.226.
- Mean jitter ARI on Markov/Ward at n=480 is .474 / .971 / .971 / .570. On the imbalanced panel it is .479 / .823 / .823 / .730. On symmetric/average at p=72, A is more stable (.844) than B (.634), C (.692) or rule-out (.408). Identical single-cluster outputs give jitter ARI 1 at gap 2.5 despite zero truth ARI.
- The original Gaussian single-population control splits in 3/4 symmetric/average and 4/4 Markov/Ward A fits. Rule-out splits in 1/4 and 3/4. B/C split in 0/4 on either tree. These finite controls demonstrate the observed behavior. Four independent seeds per regime provide little precision for rates, especially when no events occur.
- The median full fit runtime is 0.484–0.489 s across arms, including geometry, branch fitting and traversal, with one BLAS thread. The archived CHOIR fit runtimes came from different thread settings and exclude diffusion-distance construction.

## Reproducibility

Protocol, executable scripts, per-fit summaries, paired differences, per-group recovery, cell assignments, relationship strata, and SHA-256 manifests are in `results/analyses/root_validation_followup_20261109/metrics/`. One focused integration test passed. All 448 simulation rows and 224 stability pairs were independently checked against saved assignments and regenerated truth. All 117 shared live audit source/input hashes matched, and the recorded source hashes stayed fixed through verification. The four dispatcher file checks use their additional mid-run hash snapshot. No count-depth preprocessing experiment is included in this metrics panel.

## Frozen protocol

Frozen before reading follow-up results. Descriptive finite panel, no threshold tuning on CHOIR or on simulations.

CHOIR uses all 48,879 cells from the saved metadata, six saved A/B/C assignments and the two saved production rule-out assignments from the ceiling-3 rerun. Check exact identity, order, Cell_line labels, count, ARI, NMI, legacy 90/90, published CHOIR_1 ARI, and equality with both earlier saved rule-out assignments. Record input SHA-256 and compare live audit manifest entries. The saved runtime excludes precomputing diffusion distances and came from different thread counts across archives.

Use truth rows and predicted-cluster columns in the contingency table. Match one-to-one by maximizing total F1, with zero scores for unmatched truth rows. Report per-group precision, recall and F1, equal-weight macro mean, rare-group macro mean (CHOIR smallest ceil(K/4) with lexical ties, synthetic explicitly generated 5% group), and matched 90/90. Separately reproduce the legacy largest-overlap/column purity 90/90 convention. Report ARI, AMI (arithmetic normalization), homogeneity, completeness, signed cluster-count error, and pair splitting/merging, the shares of same-truth pairs assigned apart and same-predicted pairs assigned different truths, respectively. A zero denominator is reported as missing. These measure truth agreement, not the probability of a false statistical rejection.

Fresh panel uses independent RNG seeds 731901–731904, shared across all arms for each regime. Base regime uses 240 cells in four equally sized Gaussian groups, 12 variables, group centers on the vertices of a square with neighbor gap 5, independent unit-variance noise. Change one factor at a time: n=480; p=72; proportions 55/25/15/5; gap=2.5; anisotropy with covariance diag(4,1/4,1,...) rotated by 45 degrees in the two signal coordinates. Add an independent-unit-Gaussian single-population control at base n/p. Group labels are shuffled together with rows. Stability uses a fixed independent per-cell feature perturbation N(0,0.1^2), rebuilt through the entire pipeline. Compare fit quality across fresh seeds and ARI between the original and perturbed partition for the same rows.

All arms use registered symmetric/average or Markov/Ward geometry and branch fitting, alpha 0.01 and edge alpha 0.001, and run through the canonical dispatcher. A and rule-out use current production code. B/C swap only the audit hook and inherit A's empirical eligibility; C reuses B's selected dimensions and rho. One process, BLAS1. No data split. Report every fit and regime, paired deltas relative to A, and runtime including geometry and traversal. Rank methods descriptively on macro F1, AMI and ARI, stability and signed count errors; report disagreements and tradeoffs rather than declaring a winner from a single metric. The null-control sample size does not support a calibration claim.

Saved split relationships: compare B and C raw scores at alpha 0.01 by parent size, smaller-child share and k. Join same tree and node to prior spectral and distance traces only after checking size and smaller-child size. Measure multiplicity and centroid-distance tie fractions descriptively in that intersecting subset. Derive topological depth from saved linkage under the same check. No inference using node independence.
