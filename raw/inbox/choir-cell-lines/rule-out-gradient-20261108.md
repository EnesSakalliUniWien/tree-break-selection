# The gradient behind the split-count rule-out

Analysis on 2026-11-08. Local folder: `results/analyses/rule_out_gradient_20261108` (git-ignored); protocol,
deviation, code and outputs in `rule-out-gradient-scripts-20261108.py.txt`.

## Observed gradient (continuous benchmark, `rule-out-continuous-benchmark-20261108.md`)

Ordered by the cells of a two-group subtree divided by the features: 0.001 (1 of 4 groups recovered), 0.44 (1 of
6), 0.75 (3 of 4), 1.0 (2 of 3), 1.1 (4 of 4), 1.6 and above (all). Every failure lies where the subtree to be
split has at most p + 1 cells.

## Grid (frozen protocol)

Two balanced Gaussian groups, unit noise, separation delta along one feature; p in {10, 30, 100}, n from 20 to
640, delta from 0 to 6, 40 replicates (8,640 runs); the true split is tested with production
`split_count_rule_out`; count pass = union bound over 2^(n-1) - 1 splits at 0.01 with known noise.

## Three relationships

1. Invariant direction (spiked-covariance transition, Baik, Ben Arous and Peche 2005, Ann. Probab. 33:1643-1697;
   threshold 1 + sqrt(p/n) from secondary sources): share of runs with an invariant direction 21% (z in -2 to -1),
   28% (-0.5 to 0), 53% (0 to 0.5), 87% (0.5 to 1), 99% (above 1), with z = log(pi(1 - pi) delta^2) -
   0.5 log(p/n). Gets easier with more cells, as predicted.
2. Count of all splits: pass 0% at pi(1 - pi) delta^2 = 1.0, 20% at 1.56, 66% at 2.25, 95% at 4, 100% at 9; the
   prediction was about 1.39, independent of n.
3. Separability (a flat cut with every cell on its side): binding, and it gets harder with more cells. Observed
   separable share against Phi(delta/2)^n (every cell on its own side of the midpoint), p = 10, delta = 4: n 20
   0.63 vs 0.63, 40 0.58 vs 0.40, 60 0.40 vs 0.25, 80 0.23 vs 0.16, 120 0.10 vs 0.06, 160 0.05 vs 0.03; delta = 6,
   n = 640: 0.60 vs 0.42 (p = 10), 0.73 (p = 30), 0.68 (p = 100). The cut can use other directions, so observed
   shares sit somewhat above the prediction. The equidistant-cell rule is not binding here (median tie share
   about 0 among separable runs).

Plus the implementation edge: no invariant directions when n <= p + 1, so kept = 0 there (E4, all 1,600 runs),
which the theory does not require. Noise (delta = 0, random half splits): kept 0 of 1,080.

Logistic fit of kept (n > p + 1, delta > 0, 6,160 runs): with z alone pseudo-R2 0.43; adding log(n/p) gives
pseudo-R2 0.74 with a negative coefficient (-3.34, SE 0.14): beyond the transition variable, more cells per
feature lower the chance of being kept.

## Reading

The rule-out accepts a split only if (1) its direction rises above the noise band, pi(1 - pi) delta^2 > about
sqrt(p/n), easier with more cells; (2) it beats the count of all splits, pi(1 - pi) delta^2 > about 2, the same
for any n; and (3) the two sides do not overlap at all, which needs the expected number of cells on the wrong
side, n Phi(-delta/2), to stay below about 1. The third term makes the required separation grow with the number
of cells: delta about 2 Phi^-1(1 - 1/n), that is 4.7 noise SDs at 100 cells, 6.2 at 1,000 and 8.2 at 48,879. It
is also what removes noise splits (one population always overlaps). So the gradient is not only cells per
feature: the rule asks for zero overlap rather than a significant separation, and that demand rises with
sample size. In the single-cell data, the lines it recovers are lines that do not overlap their neighbours at
their subtree's size; groups that overlap even slightly are refused as the subtree grows.
