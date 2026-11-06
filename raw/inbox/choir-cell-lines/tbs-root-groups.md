# Groups directly attached to the root in the four cancer cell-line trees

On 2026-11-04, we expanded the final merge of each saved linkage to its leaf
members and joined those rows to the frozen cell metadata and diffusion graph.
All four pairs of groups partition the same 48,879 cells. The smaller groups
are pairwise disjoint. Each smaller group contains a single annotated cell
line, and each larger group contains at least 189 lines.

| Tree | Smaller child | Cells from the line | Larger child size / lines | `CHOIR_1` counts in smaller child | Diffusion graph edges across split |
| --- | --- | ---: | ---: | --- | ---: |
| Symmetric / average | `L43204`, ASPC1_PANCREAS | 1 / 290 | 48,878 / 190 | 18 (1) | 88 |
| Symmetric / Ward | `N97689`, A204_SOFT_TISSUE | 250 / 250 | 48,629 / 189 | 148 (130), 170 (120) | 0 |
| Markov / average | `N97677`, OAW28_OVARY | 151 / 161 | 48,728 / 190 | 103 (151) | 2 |
| Markov / Ward | `N97606`, NCIH727_LUNG | 267 / 270 | 48,612 / 190 | 80 (170), 217 (97) | 9 |

- The ASPC1 singleton has 88 graph neighbors, all from ASPC1 and all on the
  larger side. The recorded baseline failed both root edge gates. Its adjusted
  p-values were 0.9167 and 0.7246, respectively.
- The A204 group is the entire 250-cell disconnected graph component. Both of
  its published CHOIR clusters are wholly inside this group. The other small
  disconnected component, the 185 IM95_STOMACH cells, stays on the larger side
  in all four trees.
- The OAW28 group is exactly published CHOIR cluster 103. The 10 OAW28 cells
  on the larger side belong to clusters 86 (1) and 249 (9). One of its two
  boundary graph edges connects two OAW28 cells.
- The NCIH727 group leaves three cells of that line on the larger side, all
  from cluster 80. All nine boundary graph edges connect NCIH727 cells.
  Cluster 217 also includes one cell from another line on the larger side.

## PCA and NNLS at the root

Distances use the NNLS target, mean squared standardized Euclidean distance
across 30 PCs. We computed the all-pairs means from the group means and
variances. The cross-pair mean is the mean squared difference between the
two centroids plus the two mean population variances. The within-group mean
uses distinct pairs. Graph edge counts count each crossing edge once.

| Tree | Mean target within smaller child | Mean target across root | Mean target within larger child | Sampled pairs across root | Root length sum, native / independent | Cross-pair fit RMSE, native / independent |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Symmetric / average | — | 2.259 | 2.000 | 1 | 0.270 / 0.160 | 0.0414 / 2.33e-9 |
| Symmetric / Ward | 0.432 | 7.034 | 1.948 | 567 | 2.981 / 2.859 | 0.302 / 0.144 |
| Markov / average | 0.420 | 2.299 | 1.998 | 326 | 0.002 / 0.087 | 0.271 / 0.157 |
| Markov / Ward | 0.648 | 6.029 | 1.956 | 574 | 1.402 / 1.187 | 0.312 / 0.192 |

The two Ward groups have the largest separation in the original PCA target.
OAW28 has only two diffusion graph boundary edges, while its average PCA
distance across the root is close to the average within the larger child.
The graph boundary and PCA distance therefore describe different aspects of
the separation.

For every topology, each of the two root-edge columns equals the indicator
for a sampled pair crossing the root. The predicted distances depend on the
sum of those two lengths. The linkage ratio determines how the fitted sum is
divided between the edges. The baseline independent fit nearly interpolates
its one sampled crossing pair. Its target value, 1.424, is below the mean
2.259 across all 48,878 crossing pairs. The original fit has a recorded TBS
gate result only for the baseline. The alternative native values are the saved
iteration-limit iterates, and the independent lengths are diagnostic fits.

## Evidence

- `tbs-root-groups-evidence.json` contains the eight group summaries, the
  intersection matrix and the NNLS measurements.
- `tbs-root-groups-composition.csv` contains complete line, cancer type, CHOIR,
  source library, pool and graph component counts for both children.
- `tbs-root-groups-script.py.txt` records the calculation. It checks saved
  input hashes, the canonical cell order, exact partition coverage, graph
  components, both root-edge columns, reconstructed targets and the previously
  applied baseline root lengths.
- The local run is `results/analyses/choir_root_groups_20261104_102259/`. It
  also retains the individual cell IDs of the four smaller groups.
