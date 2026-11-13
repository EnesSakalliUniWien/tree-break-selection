# Diffusion operator against linkage on the cancer_pc cases, today's test

Run on 2026-11-12, Europe/Vienna. Local folder: `results/analyses/cancer_pc_operator_linkage_20261112`
(git-ignored); protocol, deviations, both scripts and the complete summary output are in
`cancer-pc-operator-linkage-scripts-20261112.py.txt`. Follow-up of `cancer-pc-over-splits-20261112.md`.
Working tree on top of `bc044e07`; registered cell of `tbs_family_diffusion` with two overrides, no opt-in.
No production code was changed. Environment: Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3, uv 0.9.5.

## Question

`tbs_family_diffusion` (operator "symmetric", linkage "average") and `tbs_family_diffusion_markov_ward`
(operator "markov", linkage "ward") differ in two parameters. Which of the two makes the accepted extra splits
of today's test?

## Setup

Method `tbs_family_diffusion` with `diffusion_coordinate_operator` in {"symmetric", "markov"} x
`tree_linkage_method` in {"average", "ward"}: 4 cells x 4 cases = 16 runs, one record per internal node (7,920
records), same records and the same classification (pass-through, inside, between) as the first inspection.
Self-check passed: symmetric/average and markov/ward reproduce found clusters and ARI of the two registered
methods in 8 of 8, so the method id carries nothing besides the two parameters. The protocol was written before
the run with the hypothesis "the linkage decides".

## Results

Found clusters and ARI (6 true groups):

| Case | symmetric/average | symmetric/ward | markov/average | markov/ward |
| --- | --- | --- | --- | --- |
| hub_disconnected | 6, 0.9982 | 6, 0.9982 | 6, 0.9982 | 10, 0.4459 |
| hub_connected | 5, 0.5934 | 12, 0.2761 | 8, 0.6718 | 13, 0.2651 |
| nohub_disconnected | 6, 0.9982 | 6, 0.9982 | 6, 0.9982 | 16, 0.2677 |
| nohub_connected | 6, 0.6974 | 9, 0.4290 | 15, 0.2656 | 15, 0.2344 |
| Mean ARI | 0.8218 | 0.6754 | 0.7334 | 0.3033 |

Accepted splits over the 4 cases:

| | symmetric/average | symmetric/ward | markov/average | markov/ward |
| --- | ---: | ---: | ---: | ---: |
| between groups | 16 | 16 | 16 | 17 |
| inside, all cells from one true group | 0 | 0 | 0 | 13 |
| inside, cells from several true groups | 3 | 13 | 14 | 20 |
| pass-through nodes | 3 | 1 | 1 | 0 |

All inside splits with cells from several true groups are in the two connected cases. The reference refuses
(p > 0.01) every accepted inside split in all 4 cells (3, 13, 14, 33).

Nodes with at least 40 cells, all from one true group:

| | symmetric/average | symmetric/ward | markov/average | markov/ward |
| --- | ---: | ---: | ---: | ---: |
| Nodes | 48 | 29 | 45 | 24 |
| Share of the smaller side, median | 0.124 | 0.228 | 0.206 | 0.382 |
| Smaller side under 10 cells | 21 | 3 | 14 | 1 |
| `rho` in 30 features, median | 1.214 | 1.369 | 2.947 | 3.872 |
| Edge significant for at least one child | 0 | 0 | 2 | 17 |
| Edge significant for both children | 0 | 0 | 0 | 17 |
| Sibling p-value present | 0 | 0 | 2 | 17 |
| `Sibling_BH_Different` | 0 | 0 | 0 | 12 |
| Reference p <= 0.01 | 0 | 0 | 0 | 0 |

The same 64 cells of group 4 in the disconnected cases: symmetric/average 10|54, `rho` in 30 features 0.334,
edge p (BH) 0.99; symmetric/ward 10|54, 0.334, 0.99; markov/average 12|52, 2.95, edge p 0.0041 and 7.4e-05
(one child significant), sibling p 0.55 and 0.48; markov/ward 29|35, 4.51, edge p 0.00021 for both children,
sibling p 0.026 (with hub) and 0.0012 (without).

## Reading

- The hypothesis "the linkage decides" is not supported. Splits inside one pure true group are accepted only
  when both parameters are changed (0, 0, 0, 13), and the disconnected cases over-split only in markov/ward.
- Each parameter alone already adds accepted splits inside the mixture of groups 0, 4 and 5 of the connected
  cases: 3 with neither, 13 with ward alone, 14 with markov alone, 20 with both.
- Each parameter alone makes the splits inside one true group more balanced (0.124 to 0.228 and 0.206); the
  operator alone raises `rho` in 30 features (1.214 to 2.947), the linkage alone does not (1.369).
- The edge gate inside one true group opens for both children only with both parameters (17 of 24 nodes);
  with markov alone it opens for one child at 2 of 45 nodes, and the sibling test then does not reject.

## Deviation

The self-check "split nodes equal found clusters minus 1 plus pass-through nodes" stopped the first run at the
15th fit. The identity is not general; it was weakened to an inequality and all 16 fits were rerun. The per-cell
table (Part 3 of the summary) was added after the first output.

## Not checked

Why the edge test is significant inside one Gaussian group with markov/ward (branch lengths, selection of the
split) was not examined. One seed per case, 4 synthetic cases, nodes with at least 40 cells; `rho` in 30
features is not comparable between nodes of different size. No real data.
