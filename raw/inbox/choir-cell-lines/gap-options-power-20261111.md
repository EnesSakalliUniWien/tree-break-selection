# Power of four gap options on the layouts the Gaussian best-flat-cut reference refuses

Run on 2026-11-11, Europe/Vienna. Local folder: `results/analyses/gap_options_power_20261111` (git-ignored);
protocol (written before the run), script and the summary are in `gap-options-power-scripts-20261111.py.txt`.
Simulation only; no production code was changed and nothing was promoted.

## Question

The benchmark run of the same day (`best-flat-cut-reference-benchmark-20261111.md`) showed that `rho` against
the Gaussian best-flat-cut reference refuses a group of 64 cells at offset 3 next to 256 cells, and a group at
offset 20 whose complement holds a mirrored group. Which of four additions to `rho` accepts these splits, and
what does each cost on one non-Gaussian population?

## Setup

One direction (k = 1, as in the refused cancer_pc subtree), unit within-group noise, seed 20261111, Python
3.12.11, NumPy 2.3.4. Statistics, all unchanged by shift and scale of the axis:

- rho: `pi(1-pi) D^2` of a split (the current statistic).
- mixture: twice the log-likelihood gain of a two-Gaussian mixture with common variance (EM, 5 starts) over one
  Gaussian. A statistic of the sample, not of a split.
- dip: distance of the empirical distribution from the closest unimodal one, own implementation from the
  definition, not checked against a published implementation; self-checks: 1/(2n) on an evenly spaced sample,
  0.25 on two tight distant clusters. A statistic of the sample.
- valley: log((min(c_left, c_right) + 0.5) / (c_mid + 0.5)), cell counts in three bands of width one third of
  the distance between the side means, centred at the side means and half-way between them.
- shape: minus the size-weighted mean of the two sides' skewness toward the other side.

Reference: 2,000 one-Gaussian samples per n; for rho, valley and shape the largest value over all threshold
cuts with at least 5 cells per side; p = (1 + #{reference >= value}) / 2,001. 99% points at n = 320 / 384:
rho 2.171 / 2.089, mixture 10.22 / 10.01, dip 0.0202 / 0.0179, valley 0.110 / 0.071, shape 0.130 / 0.128.

## Results: samples with p <= 0.01 / p <= 0.05

True-group split, of 1,000 samples per layout:

| Layout (cells at position) | n | rho | valley | shape |
| --- | ---: | --- | --- | --- |
| A: 256 at 0, 64 at 3 | 320 | 0 / 0 | 4 / 8 | 158 / 350 |
| B: 256 at 0, 64 at -3, 64 at 3; group at 3 against rest | 384 | 0 / 0 | 0 / 0 | 0 / 0 |
| C: 256 at 0, 64 at -20, 64 at 20; group at 20 against rest | 384 | 0 / 0 | 1,000 / 1,000 | 0 / 0 |
| D: 160 at 0, 160 at 3 | 320 | 643 / 839 | 855 / 947 | 183 / 358 |
| E: 256 at 0, 64 at 4 | 320 | 932 / 986 | 546 / 698 | 154 / 316 |

Whole sample (mixture, dip) and best threshold cut (rho, valley, shape), of 1,000 samples per layout:

| Layout | best rho | mixture | dip | best valley | best shape |
| --- | --- | --- | --- | --- | --- |
| A | 278 / 603 | 1,000 / 1,000 | 3 / 39 | 27 / 68 | 1,000 / 1,000 |
| B | 0 / 0 | 3 / 29 | 5 / 20 | 0 / 5 | 9 / 59 |
| C | 0 / 0 | 1,000 / 1,000 | 1,000 / 1,000 | 1,000 / 1,000 | 0 / 0 |
| D | 1,000 / 1,000 | 1,000 / 1,000 | 664 / 860 | 979 / 999 | 1 / 44 |
| E | 1,000 / 1,000 | 1,000 / 1,000 | 247 / 455 | 888 / 953 | 1,000 / 1,000 |

One population, n = 320, of 2,000 samples (20 / 100 expected for the Gaussian):

| Population | best rho | mixture | dip | best valley | best shape |
| --- | --- | --- | --- | --- | --- |
| Gaussian, fresh samples (level check) | 6 / 71 | 10 / 72 | 15 / 102 | 13 / 84 | 22 / 124 |
| Gamma, shape 4 | 65 / 249 | 2,000 / 2,000 | 15 / 107 | 15 / 83 | 2,000 / 2,000 |
| Student t, 5 degrees of freedom | 0 / 0 | 1,403 / 1,700 | 4 / 47 | 0 / 3 | 30 / 63 |

The median cut of the same one-population samples is rejected in at most 3 / 49 samples (rho, Gaussian) and in
0 to 2 samples for all other cells.

## Reading

- B, the connected cancer_pc layout, is not accepted by any option: at most 9 of 1,000 samples at 0.01. Along
  its axis the layout has valleys at 0.999 of the neighbouring mode, skewness 0 and kurtosis 3.000 (density
  calculation in the benchmark capture), so it is close to one Gaussian.
- A (64 cells at offset 3 next to 256) is accepted only by the mixture (1,000 of 1,000) and by the shape of
  the best cut (1,000). Both also reject every gamma sample, and the mixture rejects 1,403 of 2,000 Student t
  samples: they detect that the sample is not Gaussian, not that it holds two groups.
- C (mirrored distant group) is accepted by dip, valley and mixture in 1,000 of 1,000 samples and by rho and
  shape in none.
- Dip and valley keep the level on the gamma and Student t populations (0 to 15 of 2,000 at 0.01).
- On D the valley of the true-group split has more power than rho (855 against 643); on E it has less (546
  against 932). The dip has 664 and 247.
- Against the expectation written before the run: confirmed for dip and valley on A, B and C, for rho, and for
  B; not confirmed for D, where dip and valley were expected near 1,000 (664 and 855), and for A, where shape
  also exceeds the level.

## Consequence for the benchmark panel

No option would change the benchmark result. The connected cancer_pc cases correspond to B, which nothing
accepts. The disconnected cases correspond to C, where the pass-through traversal already gives the correct
partition (ARI 0.998).

## Limits

- One direction only; k of 2 and above, the tree's choice of split and complete fits were not simulated.
- For valley and shape the largest value over threshold cuts is not proved to bound every split.
- The valley's band width (one third) is a chosen value; other widths were not tried.
- The dip implementation was not compared with a published one; because its reference is simulated with the
  same code, the comparison of power does not depend on that.
- The reference has 2,000 samples; the level check shows its noise (6 to 22 at 0.01 where 20 are expected).
- A rule that accepts a split when either rho or a gap measure passes needs its level divided between the
  two; that rule was not simulated.
