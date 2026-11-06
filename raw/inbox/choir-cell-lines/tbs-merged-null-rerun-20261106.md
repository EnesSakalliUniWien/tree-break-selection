# Count-null rerun after the singleton-peel merge

Run on 2026-11-06 on merge commit `50455316`, after
`tbs-merged-rerun-20261106.md`.

## Protocol

The 600 homogeneous count nulls of the frozen 2026-11-03 recovery protocol were
rerun with all four variants: 300 with library log-sigma 0.25 and 300 with 1.0,
400 cells, 3,000 negative-binomial genes. `validate.py`, `compare.py` and
`protocol.json` were copied unchanged from
`results/analyses/choir_tbs_recovery_comparison_20261103_174622/`;
`validate.py` matches `tbs-recovery-validate-script.py.txt` and the recorded
SHA-256 `2f0448b2...`. Six shards ran in parallel (model by replicate modulo
three), as in the original execution. The self-check returned one cluster.

All 2,400 fits completed in 570 s with no failures or tracebacks. Every count
dataset's SHA-256 matches the original study, dataset by dataset. The endpoint
and bounds follow the original `aggregate.py`: a false partition is final
K > 1, with a one-sided 95% exact binomial upper bound, failures counted as
errors, against the frozen target 0.01.

## Results

False final partitions out of 300 nulls, before and after the merge:

| Variant | Low depth variation | High depth variation |
| --- | --- | --- |
| Symmetric / average (production baseline) | 0 → 2 | 262 → 264 |
| Symmetric / Ward | 295 → 295 | 228 → 228 |
| Markov / average | 34 → 41 | 298 → 298 |
| Markov / Ward | 297 → 297 | 300 → 300 |

The production baseline's low-depth upper bound rises from 0.00994 to 0.020836,
so it no longer meets the 0.01 target. No variant meets the target under
both null models. In paired comparisons on the baseline, two datasets per null
model newly split and none returned to one cluster; median K is unchanged
(1 low, 5 high).

Both newly split low-depth nulls (replicates 90 and 178) have a 1/399 singleton
root; traversal passes through two or three singleton peels, then makes one
split, giving K = 3 and 4. 108 of the 300 low-depth baseline nulls have a
singleton root, and only these two used a pass-through. The two newly split
high-depth nulls (138 and 246) have 74/326 and 76/324 roots, with one
pass-through each.

## Interpretation limits

Two of 300 is a point rate of 0.67%. With 300 nulls any single false partition
puts the upper bound above 0.01, so the result shows the target is no longer
established, not that the rate clearly exceeds 1%. Signal mixtures were not
rerun.

## Files

- `tbs-merged-null-summary-20261106.csv`: per-model, per-variant summary.
- Local-only outputs: `results/analyses/choir_count_null_rerun_merged_20261106/`.

SHA-256:

- `validation_null_summary.csv` 0a51de3525c35c91c5496b3accd530c53aa24243a8aca11f6411c552314a0dc8
- `validation_dataset_results.csv` 76a1c8bfef3f65c935bb8818fdd307c76998ea04d0027612bea95001e428de69
- `validate.py` 2f0448b298a64caae61a9f48e73c49d45201bf293800b77af23c5c2232d27017
- `compare.py` 05d998d09673480c5346167eaeb3bf068388148705fad63f4232fc3ad268f4cb
- `protocol.json` 0378e282b79fd2ee23c6ca693f66e760cedc48eb22d4cb8a4ee074db085d4044
