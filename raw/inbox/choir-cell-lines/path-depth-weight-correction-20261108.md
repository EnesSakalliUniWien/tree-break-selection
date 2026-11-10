# Large-split sibling correction with tree depth, the weight and whole paths

Analysis on 2026-11-08, code at `cbc317b2` plus the staged working tree. Local folder:
`results/analyses/path_depth_weight_correction_20261108` (git-ignored); protocol, deviations and scripts
captured in `path-depth-weight-correction-scripts-20261108.py.txt`. Protocol frozen before the full run;
one descriptive endpoint added after a 2-dataset pilot (see Deviations).

## Question and definitions

Does the sibling correction come closer to what no-group data need for large splits, without losing real
groups, when it takes into account tree depth, the reference weight and whole paths? User definitions: tree
depth = number of splits between the root and the split; whole paths = every root-to-leaf path through the
split (its ancestors and its descendants).

## Candidates

All use only the tree's own reference records and today's estimator (weighted ratio, Gaussian kernel,
weighted-SD bandwidth, family-baseline cap, c >= 1, p = chi2.sf(T/(s*c), df)). Production code runs end to
end (edge gate, Tree-BH, traversal BH, pass-through); only the factor c of each tested split is swapped.

| Name | Change from today |
| --- | --- |
| today | none; replay equals production c on 11,488/11,488 tested splits |
| depth | tree depth as a third kernel axis |
| weight | only edge-blocked reference records (kept by a closed gate above them, weight 1) |
| paths | only reference records that are ancestors or descendants of the split |
| all | the three together |

Data: the 900 one-group count datasets of 2026-11-07 (400 cells; log-depth sd 0.25 `low`, 0.42
`matched_median`, 1.0 `high`; 300 each) and 35 datasets with real groups (`signal_balanced`,
`signal_imbalanced`, `signal_singleton`: 1,000 cells, 4 groups, 10 each; `signal_density_variation`: 3,000
cells, 4 groups, 5). Both trees, production alphas (edge 1e-3, sibling 1e-2).

## E1 Tested splits of 40 or more cells on no-group data

Share with adjusted p <= 0.01 (target about 1%), median applied factor in brackets.

| Tree | sd | today | depth | weight | paths | all |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Markov/Ward | 0.25 | 71.8% (2.5) | 71.9% (2.6) | 66.4% (2.7) | 66.7% (2.8) | 49.9% (3.4) |
| Markov/Ward | 0.42 | 73.1% (2.5) | 72.7% (2.5) | 67.6% (2.7) | 67.7% (2.7) | 51.2% (3.3) |
| Markov/Ward | 1.0 | 84.8% (2.3) | 85.3% (2.2) | 81.3% (2.4) | 76.4% (2.5) | 59.6% (3.3) |
| Default | 0.42 (163 tested) | 54.6% (2.4) | 39.9% (2.7) | 53.4% (2.5) | 50.3% (2.5) | 37.4% (2.9) |
| Default | 1.0 | 91.6% (1.9) | 89.2% (1.9) | 88.6% (2.1) | 86.2% (2.0) | 80.0% (2.4) |

By size, Markov/Ward sd 0.25: 40-159 cells needed 7.0, applied 2.5 (today) to 6.1 (all), pass 61% to 33%;
160-400 cells needed 7.6, applied 2.6 to 2.9, pass 94% to 85%. 10-39 cells: today 1.0% pass; paths and all
raise it to 4.4% and 5.4%.

## E2 False partitions on no-group data (of 300)

| Tree | sd | today | depth | weight | paths | all |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Markov/Ward | 0.25 | 298 | 295 | 290 | 300 | 292 |
| Markov/Ward | 0.42 | 296 | 291 | 291 | 299 | 291 |
| Markov/Ward | 1.0 | 300 | 300 | 298 | 300 | 299 |
| Default | 0.25 | 2 | 4 | 2 | 1 | 3 |
| Default | 0.42 | 50 | 37 | 50 | 51 | 37 |
| Default | 1.0 | 259 | 260 | 258 | 260 | 259 |

Median clusters on Markov/Ward: 8 (today) to 6 (all).

## E3 Data with real groups (mean ARI)

| Tree | Signal | today | depth | weight | paths | all |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Default | balanced | 0.800 | 0.800 | 0.800 | 0.800 | 0.800 |
| Default | imbalanced | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| Default | singleton | 0.880 | 0.880 | 0.880 | 0.880 | 0.880 |
| Default | density variation | 0.651 | 0.651 | 0.651 | 0.650 | 0.651 |
| Markov/Ward | balanced | 0.315 | 0.312 | 0.319 | 0.325 | 0.445 |
| Markov/Ward | imbalanced | 0.079 | 0.080 | 0.090 | 0.096 | 0.220 |
| Markov/Ward | singleton | 0.092 | 0.092 | 0.093 | 0.119 | 0.210 |
| Markov/Ward | density variation | 0.025 | 0.025 | 0.026 | 0.028 | 0.044 |

Median clusters on Markov/Ward balanced: 20 (today) to 14.5 (all); true 4.

## Decision rule (fixed before results)

- Fixes large splits (Markov/Ward, sd 0.25 and 0.42, tested splits >= 40 cells pass <= 2%): no candidate.
- Keeps real groups (default tree, mean ARI within 0.02 of today on every signal model): all candidates.

## E4 Why: reference scores against tested scores at the same size

Score = T/(s*df), about 1 for one population without selection. Reference score = sum(w*T)/sum(w*s*df);
tested score = median over tested splits. Markov/Ward:

| sd | Size | Own edges tested (n, score) | Edge blocked (n, score) | Tested splits |
| --- | --- | ---: | ---: | ---: |
| 0.25 | 10-39 | 1,702, 1.9 | 12,806, 7.9 | 14.0 |
| 0.25 | 40-159 | 1,299, 2.3 | 442, 3.2 | 4.8 |
| 0.42 | 10-39 | 1,598, 2.0 | 12,974, 8.2 | 13.5 |
| 0.42 | 40-159 | 1,312, 2.2 | 686, 3.2 | 5.1 |
| 0.42 | 160-400 | 11, 3.0 | 8, 5.4 | 7.8 |

Reference splits of 2-9 cells score 12 to 15 (edge blocked). Default tree, sd 0.42: 40-159 references score
1.3 and 1.8, tested 3.9; 160-400 references 1.3 and 3.1, tested 4.2.

At every size the tested splits score about 1.3 to 7 times more than the reference splits, although every
split here is a no-difference split. The edge gate that decides "tested" or "reference" compares each side
with the whole, which measures nearly the same difference as the sibling score, so the gate sorts splits by
their score. Depth, weight and path choose among reference splits, but every reference split is on the
low-score side of that sort, so no reweighting of them reaches the needed factor for large splits. This is a
reading of the measured scores and the code; it was not tested separately.

## Deviations

- Pilot on `low` 0 and `signal_balanced` 0 before the full run: replay 40/40; the `weight` candidate moved
  the Markov/Ward large-split factor from 2.25 to 2.28 (378 of 388 records already edge blocked). E4 was
  added because of this, before the full run. No candidate, rule or threshold changed.
