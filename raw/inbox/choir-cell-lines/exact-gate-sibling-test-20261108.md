# Exact sibling test given that the edge gate opened

Analysis on 2026-11-08, code at `cbc317b2` plus the staged working tree. Local folder:
`results/analyses/exact_gate_sibling_test_20261108` (git-ignored); protocol, deviations, scripts and the
summary output are in `exact-gate-sibling-test-scripts-20261108.py.txt`. It replaces the approximation in
`gate-conditioned-sibling-test-20261108.md` (see the correction there).

## Conditioning

`split_edge_test_rows` gives each child a test half (edge contrast) and a basis half (edge PCA basis). The
test conditions on both basis halves, the test-half pooled mean, the direction of the test-half contrast and
the plug-in quantities production already holds fixed. The free quantity is the test-half length R (factor f,
f = 1 observed). Both edge statistics are f^2 times their observed value, so the Tree-BH gate is open exactly
for f >= f*. The sibling statistic T(f) is production's `sibling_divergence_test` on the moved all-cell means.
Under no difference R^2 ~ c_t chi2(p), p = 30; p-value = P(T(f) >= T_obs and f >= f*) / P(f >= f*).

Self-checks over all 935 datasets: edge statistics recomputed with production code equal f^2 times observed
at f = 0.5 and 2 (22,976/22,976); T(1) equals the production statistic (11,488/11,488); test rows equal the
edge test's rows (479,500/479,500 edges); Tree-BH copy (1,325/1,325 fits); gate open for tested splits
(11,488/11,488) and closed for own-tested references (12,296/12,296); replay (11,488/11,488);
{T(f) >= T_obs, f >= f*} = [1, infinity) for all but 1.2-1.7% of 10-39-cell splits and every larger split.

Arms: `exact_gate` (c_t = 1), `exact_gate_ct` (c_t from the references' R^2 with today's kernel, weights and
cap), `exact_gate_both` (c_t by truncated likelihood, own-tested references truncated at their exact gate).

## Results (no-group data, tested splits with p <= 0.01)

| Tree | sd | Size | today | exact_gate | exact_gate_ct | exact_gate_both |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| Markov/Ward | 0.25 | 10-39 | 1.0% | 27.9% | 12.6% | 12.6% |
| Markov/Ward | 0.25 | 40-159 | 61.0% | 69.5% | 63.7% | 63.7% |
| Markov/Ward | 0.25 | 160-400 | 94.3% | 98.8% | 98.8% | 98.8% |
| Markov/Ward | 0.25 | >= 40 | 71.8% | 79.0% | 75.1% | 75.1% |
| Markov/Ward | 0.42 | >= 40 | 73.1% | 83.3% | 77.4% | 77.4% |
| Default | 0.42 | >= 40 (163) | 54.6% | 55.8% | 55.2% | 55.2% |

False partitions of 300: Markov/Ward 298/296 (today, sd 0.25/0.42) to 300/300 (every exact arm); default tree
2 to 3 (sd 0.25) and 50 to 57 (sd 0.42). Real groups: default-tree mean ARI unchanged or higher (density
variation 0.651 to 0.663), Markov/Ward within 0.02.

Decision rule (frozen): no arm fixes large splits; all keep real groups.

## Why

Median R^2/p of tested no-group splits (the test-half length the tree's choice inflates): 2.0-2.2 for 10-39
cells, 3.1-3.7 for 40-159, 4.4-7.3 for 160-400 (Markov/Ward). The references give c_t of only 1.0-1.3: their
R^2/p is about 1 or below (0.63 for 2-9-cell references on `low` 0), and for splits larger than every
reference the production cap sets c_t to the reference baseline, which is below 1. The gate cut (median f*^2
0.27-0.64 for >= 40 cells) removes little. So with the basis half held fixed the test is about the test-half
length, which the tree's choice inflates 2 to 7 times, and nothing in the reference splits measures that. The
improvement seen with the approximation came from pairing a too-low cut with today's c, which is estimated
for the all-cell statistic.
