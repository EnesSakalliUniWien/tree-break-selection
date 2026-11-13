# Recovery by statistical profile

This is a descriptive replay of saved assignments for 48,879 CHOIR cells and four unperturbed synthetic imbalance seeds. It recomputes truth group profiles from the archived 30-PC CHOIR matrix and the regenerated 12-feature synthetic inputs. Cell IDs, row order, truth labels, and all 1,648 existing per-group recovery rows were checked. It also measures 209 explicit CHOIR subtrees selected by the intersection of the prior B/C tested-node, spectral, and distance traces. This analysis fits no clustering method.

AMI is the primary agreement score in the research plan. The AMI entries below come from the verified saved metrics panel. Per-group precision, recall and F1 use a single one-to-one assignment maximizing total F1. Unmatched truth groups score zero. The rare summary uses the 48 smallest CHOIR lines (56–145 cells) or the predefined 12-cell synthetic component. The legacy shared macro metrics retain their existing majority remapping. Pair split share measures the fraction of same-truth pairs assigned apart. Pair merge share measures the fraction of same-prediction pairs with distinct truth. An undefined denominator yields a missing value.

| Markov/Ward panel | Test | AMI | Matched macro F1 | Matched rare F1 |
|---|---|---:|---:|---:|
| CHOIR | A | .8313 | .6173 | .7007 |
| CHOIR | B | .9346 | .7996 | .7148 |
| CHOIR | C | .9330 | .7463 | .5565 |
| CHOIR | rule-out | .8839 | .7408 | .8037 |
| Synthetic imbalance | A | .6140 | .7170 | .9345 |
| Synthetic imbalance | B | .6120 | .5095 | .0000 |
| Synthetic imbalance | C | .6120 | .5095 | .0000 |
| Synthetic imbalance | rule-out | .7077 | .7799 | .7300 |

## Measured profiles

- **Synthetic imbalance.** The 12-cell group has mean effective rank 6.595 and mean Euclidean separation 5.108 from its nearest centroid. The other groups have mean effective rank 10.176–11.472 and mean nearest separation 4.889–5.031. Thus the smallest group is not the least separated by this centroid measure, yet it is unmatched under B/C on all four Markov/Ward seeds. Its sample covariance rank is 11 in 12 features, which is the maximum possible with 12 observations. The effective-rank and covariance contrasts include finite-sample effects.
- **CHOIR truth groups.** Sample covariance is full rank in all 190 groups on the 30 archived PCs, with effective ranks from 5.714 to 20.784. On Markov/Ward, B's mean matched F1 is .484 in the lowest nearest-neighbor Mahalanobis separation quartile and .969 in the highest. Rule-out gives .546 and .888. B's matched F1 is .947 in the lowest centroid-overlap quartile and .510 in the highest. B's mean F1 also rises from .696 to .923 across the lowest and highest effective-rank quartiles. The CSV reports quartile boundaries and results for all eight fits. Profiles are stratified one at a time, so the differences do not isolate causal effects of separation, covariance, rank, or size.
- **Explicit CHOIR subtrees.** Split membership comes from the saved linkage and is checked against node and smaller-child sizes. The supplied geometry disagreement is the saved fraction of rows whose centroid-margin signs differ between tree space and test space. On Markov/Ward, its median is .005 among 75 nodes passing both B and C raw .01 thresholds and .036 among 91 passing neither. The corresponding median pooled Mahalanobis distances are 7.333 and 4.071. The one C-only node is preserved. These are nested selected nodes, and the geometry disagreement is a distance diagnostic. No population topology or statistical null is inferred from it.

The pair and subtree effective sample size is (n_L n_R/(n_L+n_R)). Covariance effective rank is the exponential entropy of the nonnegative eigenspectrum. Mahalanobis separation uses pooled within-side covariance and is missing for numerical rank deficiency. The in-sample centroid-overlap proxy averages the opposite-centroid share across the two sides and gives ties half weight. Geometry disagreement is missing unless supplied. All 209 archived subgroup contrasts have full pooled covariance rank in this panel.

## Reproduce

Run `.venv/bin/python reports/subtree_clustering_research_20261109/profiles/run.py` from the project root with one BLAS thread. The manifest hashes the saved inputs and source. The group, pair, stratum, subtree, and per-arm recovery CSVs include both registered methods and all four saved arms. These descriptive profiles do not time candidate fits or assess the research runtime cap.
