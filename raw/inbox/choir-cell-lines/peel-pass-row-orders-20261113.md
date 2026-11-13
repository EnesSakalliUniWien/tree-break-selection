# The peel pass on the Markov/Ward tree over five row orders

Run on 2026-11-13, Europe/Vienna. Local folder: `results/analyses/peel_pass_row_orders_20261113` (git-ignored);
protocol, script and outputs are in `peel-pass-row-orders-scripts-20261113.py.txt`. Follow-up of
`peel-pass-20261112.md`. No production code was changed. Environment: Python 3.12.11, NumPy 2.3.4, SciPy 1.16.3,
uv 0.9.5.

## Setup

CHOIR cells (48,879 cells, 30 PCs, 190 cell lines), Markov/Ward tree, opt-in reference on, third version of the
draft definition of k, peel pass with limit 15 at run time. The fit of 12 November is the data order of the
rows; four more fits use the row orders with seeds 1, 2, 3 and 777, the same orders as the fits without the peel
pass (`results/analyses/variant_a_fixed_20261112/choir_v3_<order>_metrics.json`). The protocol with four
predictions was written before the runs. Self-check passed in every fit: 13,641 peels of one cell and 1,867
small peels in the mapping.

## Result

| Order | Peel pass | Found clusters | ARI | NMI | Lines recovered | Split only | Merged only | Both | Largest cluster | Calls | Accepted by BH |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| data | no | 206 | 0.8359 | 0.9422 | 127 | 16 | 30 | 17 | 1,990 | 1,838 | 150 |
| data | yes | 210 | 0.8652 | 0.9442 | 127 | 17 | 29 | 17 | 1,990 | 2,210 | 153 |
| 1 | no | 211 | 0.8401 | 0.9420 | 126 | 16 | 31 | 17 | 1,990 | 1,838 | 153 |
| 1 | yes | 217 | 0.8706 | 0.9440 | 126 | 17 | 30 | 17 | 1,990 | 2,210 | 157 |
| 2 | no | 215 | 0.8411 | 0.9414 | 126 | 18 | 29 | 17 | 1,990 | 1,838 | 155 |
| 2 | yes | 217 | 0.8694 | 0.9438 | 126 | 19 | 28 | 17 | 1,990 | 2,210 | 157 |
| 3 | no | 213 | 0.8404 | 0.9415 | 126 | 18 | 29 | 17 | 1,990 | 1,838 | 153 |
| 3 | yes | 215 | 0.8687 | 0.9444 | 126 | 19 | 29 | 16 | 1,990 | 2,210 | 155 |
| 777 | no | 210 | 0.8349 | 0.9416 | 126 | 17 | 30 | 17 | 1,990 | 1,838 | 151 |
| 777 | yes | 214 | 0.8642 | 0.9437 | 126 | 18 | 29 | 17 | 1,990 | 2,210 | 154 |

Over the five orders: ARI 0.8642 to 0.8706 with the peel pass against 0.8349 to 0.8411 without; found clusters
210 to 217 against 206 to 215; lines recovered 126 or 127 in both. The gain in ARI within one order is 0.0283
to 0.0305.

## Predictions

- Z1 held: the ARI is above 0.8411 in all five orders (smallest 0.8642).
- Z2 held: in every order the ARI with the peel pass is higher than without it.
- Z3 held: found clusters 210 to 217, a difference of 7.
- Z4 held: 126 or 127 lines recovered.

## Reading

- The ARI 0.8652 of 12 November is not an effect of the row order: the two ranges do not overlap, and the gain
  within an order (0.0283 to 0.0305) is more than four times the spread between orders (0.0064 with, 0.0062
  without).
- The gain does not come from more recovered lines (126 or 127 in every fit). Per order the peel pass adds 372
  calls, 2 to 4 accepted splits and 2 to 6 found clusters.

## Not checked

- Which accepted splits carry the gain and which cell lines they separate.
- Limits other than 15 cells. The default tree in other row orders.
- Today's default sibling test with the peel pass. The edge gate itself is unchanged and not valid for a tree
  built from the tested cells. One dataset and one tree.
