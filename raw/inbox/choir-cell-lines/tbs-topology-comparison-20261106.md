# Cancer tree topologies and the Markov/Ward production method

Analysis and implementation on 2026-11-06.

## Topology comparison

`tbs-topology-comparison-script.py.txt` read the saved linkages of the merged
symmetric/average rerun and of the three 2026-11-03 alternatives (same cells,
kernel and eigenpairs; `tbs-remedy-proposal.md`). A cell line has a clean
subtree when some node holds at least 90% of the line and is at least 90% that
line. Best F1 is the line's maximum over all nodes. The 190-cluster cut uses
`fcluster(..., criterion="maxclust")` and is given the number of lines.

| Topology | Root split | Nodes >= 1,000 leaves peeling <= 10 | Lines with clean subtree | Median best F1 | ARI, 190-cluster cut |
| --- | --- | ---: | ---: | ---: | ---: |
| Symmetric / average (TBS run) | 1 / 48,878 | 0.181 of 155 | 42 | 0.558 | 0.076 |
| Symmetric / Ward | 250 / 48,629 | 0.000 of 101 | 45 | 0.564 | 0.259 |
| Markov / average | 151 / 48,728 | 0.224 of 156 | 134 | 0.975 | 0.274 |
| Markov / Ward | 267 / 48,612 | 0.019 of 106 | 134 | 0.975 | 0.731 |

TBS on the symmetric/average tree recovers 27 of its 42 clean lines. Of the
120 lines dominated by the giant cluster, 13 have a clean subtree anywhere in
that tree. Coordinates set whether lines exist as subtrees (42–45 symmetric,
134 Markov); linkage sets the peel structure. Published `CHOIR_1` labels score
ARI 0.664.

## Production change

Markov coordinates are now a `coordinate_operator="markov"` option of
`compute_diffusion_coordinates`, `adaptive_diffusion_geometry` and
`feature_family_diffusion_geometry`: `D^(-1/2)` times the symmetric
coordinates, the right eigenvectors of `D^(-1) W`. The default stays symmetric.
The opt-in registered method `tbs_family_diffusion_markov_ward` uses Markov
coordinates, Ward linkage and fixed-topology NNLS branch lengths; it is not in
`DEFAULT_METHODS`. Tests check the eigenvector relation, metadata, dispatch
forwarding and opt-in status; the full ordered suite passed 1,043 tests.

## Surrogate smoke run

On the development seed of the four `cancer_pc` surrogates (six populations):

| Case | `tbs_family_diffusion` K / ARI | Markov/Ward K / ARI |
| --- | --- | --- |
| hub, disconnected | 6 / 0.998 | 10 / 0.446 |
| hub, connected | 5 / 0.593 | 13 / 0.265 |
| no hub, disconnected | 6 / 0.998 | 17 / 0.242 |
| no hub, connected | 6 / 0.697 | 16 / 0.229 |

## Limits

In the frozen count-null validation the Markov/Ward variant split 297/300
low-depth and 300/300 high-depth nulls, and its full-data NNLS did not
converge within 1,000 iterations (`tbs-recovery-comparison.md`). The method was
not rerun on the full cancer data with the merged gates. Containing the lines
as subtrees does not make its gate decisions calibrated.
