# Native statistics across changed subtree contexts

The saved panel contains 160 successful fits on 16 reused datasets and 16,128 reference-node comparisons. Each reference tree is compared with eight alternatives, and nested nodes share observations. These are descriptive counts of comparisons. The 16 datasets are the replication unit.

## Four overlapping structural cases

| Condition | Comparison rows |
|---|---:|
| Case 1, focal clade absent | 15,316 |
| Case 2, both original child clades absent | 5,271 |
| Case 2 restricted to present focal clades | 399 |
| Case 3, focal clade present and original parent clade absent | 401 |
| Case 4, case 3 with both original child clades absent | 102 |
| At least one original child absent | 10,045 |

The 399 present-node case-2 rows include 200 roots. Case 3 asks whether the original parent clade exists anywhere in the alternative tree. The matched-context table instead asks whether the immediate parent is the same. A retained parent clade can be a higher ancestor. Counts overlap by design.

| Reference node size | Present clade | Absent clade | Total |
|---|---:|---:|---:|
| 2–3 | 270 | 8,282 | 8,552 |
| 4–7 | 11 | 3,789 | 3,800 |
| 8–15 | 14 | 2,058 | 2,072 |
| 16–31 | 202 | 790 | 992 |
| 32–63 | 59 | 397 | 456 |
| 64, root | 256 | 0 | 256 |

These rows include small terminal groupings and all other internal clades. They do not measure loss of known factor clades. The accompanying relationships analysis measures those factors directly.

## Collected sibling tests, branch multipliers and traversal

The table compares present clades only. (T) is the raw collected sibling statistic, (g) is the saved branch variance multiplier, and (Tg) is a fixed-context, unregularized scale diagnostic. Under the native ridge, multiplying the recorded statistic by (g) need not exactly recover a score recomputed at (g=1). Changes use both absolute and relative tolerances of (10^{-10}). Inflation has its own denominator of finite paired values.

| Immediate parent | Child partition | Comparisons | Changed (T) | Changed df | Changed (g) | Changed (Tg) | Changed inflation / finite pairs |
|---|---|---:|---:|---:|---:|---:|---:|
| Same | Same | 19 | 19 | 0 | 19 | 1 | 6/6 |
| Same | Changed | 84 | 84 | 1 | 81 | 82 | 8/8 |
| Changed | Same | 277 | 277 | 0 | 277 | 22 | 24/26 |
| Changed | Changed | 176 | 176 | 11 | 171 | 175 | 17/17 |
| Root | Same | 8 | 8 | 0 | 8 | 0 | 8/8 |
| Root | Changed | 248 | 248 | 0 | 248 | 248 | 243/246 |

In all 19 comparisons with the same immediate parent and children, the incoming test rows, incoming basis rows and outgoing basis rows agree. The branch factor and raw (T) still change in every pair. Across the 296 nonroot same-child comparisons, the median absolute (Tg) difference is (2.53\times10^{-11}). The 23 tolerance exceptions range from (7.68\times10^{-10}) to (3.50\times10^{-7}) in absolute units and (1.28\times10^{-10}) to (5.84\times10^{-8}) relative to the larger absolute statistic. This rescaling is a counterfactual diagnostic for the unregularized fixed-context statistic. The small residuals do not identify an independent cause beyond numerical and fitted effects. With changed children, (Tg) changes in 257/260 nonroot pairs. Those scores compare different sibling contrasts.

Among the 812 present alternative clades, 479 have an unavailable active BH decision despite a finite collected raw p-value. The remaining 331 record active rejection and two active non-rejection. Traversal visits 523 and leaves 289 unvisited. Structural absence, gate availability, visits and final boundaries are recorded separately in the JSON. A finite raw p-value in these replays does not establish a valid selected-null reference.

Two saved null examples use the same replicate 0 Markov/Ward reference.

- **N72, same children and changed parent.** The focal pair `14,31` is unchanged. Its parent changes from `14,31,39,48` to `11,14,31` in `cov_0_1`. Raw (T) changes from 4.435 to 4.464, (g) changes from 1.353 to 1.344, and df remains 1. Neither node is visited and neither has an active BH decision.
- **N104, case 4.** The focal clade `13,29,30,59` remains in `cov_5_6`. Its original parent disappears and the original children `13,30 | 29,59` are replaced by `13,29,30 | 59`. Raw (T) changes from 4.707 to 4.458 at df 2. Both gate decisions are unavailable and neither node is visited. The two statistics use different child partitions.

## Covariance invariance and inherited versus rebuilt geometry

The script independently recomputes covariance spectra for 9,018 unique dataset–clade pairs and verifies all 10,080 saved node spectra. The maximum absolute eigenvalue error is (3.38\times10^{-14}). The raw sample covariance is mathematically fixed by the same rows and feature coordinates. Guarded covariance, projected contrast, branch scaling and empirical correction can still depend on the tree.

All 32 native reference fits match the previous saved runs in dataset hash, status, root split, root rejection and partition (ARI 1). AMI agrees exactly and root p-values agree at relative tolerance (10^{-8}).

| Location | Inherited/rebuilt comparisons | Rank mismatch | Equal rank, different span | Same span |
|---|---:|---:|---:|---:|
| Root | 128 | 0 | 80 | 48 |
| Nonroot | 7,936 | 1,130 | 6,575 | 231 |

Rank mismatch can result from retaining all six inherited feature coordinates while a small child has fewer positive sample covariance directions. For `cov_full` nonroot comparisons, the median projector gap is 2.000 across all pairs and (9.33\times10^{-16}) among equal-rank pairs. For the leading two-block view, the equal-rank median gap is 1.925. A rank mismatch alone does not measure angular rotation.

No audited local covariance trace is zero. In the initial traces, the positive-rank cutoff crosses a numerical tie block in 56 comparisons and the actual rebuilt basis ends inside such a block in three. The corrected diagnostic excludes a block unless it is wholly positive under the numerical rank rule. Its saved-data sidecar recomputes all 8,064 comparisons. Only three rebuilt ranks change, from 4 to 3 for tied-population crossed replicate 0 N103 and replicate 3 N100, and from 5 to 4 for common-mean null replicate 1 N80, all under `cov_full`. No corrected rebuilt cutoff splits a numerical tie block. The counts in the table are the corrected counts and agree with the initial counts. No native fit or original trace was changed.

These comparisons concern covariance bases in common feature coordinates. Inherited parent diffusion coordinates restricted to child rows, the induced parent graph, and the refitted child diffusion graph are separate objects. In the prior saved diffusion panel, induced and rebuilt symmetric operators differed by more than (10^{-10}) in 316/320 child comparisons. The result is a paired geometry diagnostic on the same observations.

## Evidence

Run `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python -m reports.subtree_context_statistics_20261109.summarize_native` from the repository root. It reads the finished panel and writes `native_summary.json`, which contains the gate and traversal transitions, finite denominators, residual ranges, both initial and corrected rank summaries, reference parity, examples and SHA-256 values. It also writes `native_covariance_spans_corrected.jsonl.gz`. The executed native runner is preserved as `source_during_native_run.py.txt` and verified against the original native manifest. The updated helper and summarizer are hashed separately. This summary reran no clustering fits.
