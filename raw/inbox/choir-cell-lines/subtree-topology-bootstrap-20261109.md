# Fixed-leaf bootstrap and Ward tie families

Date: 2026-11-09, Europe/Vienna. The fixed protocol, runnable analysis,
source hashes and evidence are in
`results/analyses/subtree_topology_reference_20261109/bootstrap/`.
Three configurations used six explicit measured coordinates per leaf
and 1,000 replicates per scheme, for 12,000 replicates. Rooted
Robinson–Foulds distances count nontrivial clades in the symmetric
difference. The evidence also gives frequencies for every nontrivial
edge split and the histogram of the minimum rooted RF to the exact
unperturbed tie family.

| Geometry | Scheme | Distinct rooted trees | Topology entropy (bits) | Outside exact tied family | Exact comparator recovery |
| --- | --- | ---: | ---: | ---: | ---: |
| Planted six-leaf hierarchy | Feature | 6 | 0.460 | 66 | 934 |
| Planted six-leaf hierarchy | Gaussian noise, sigma=.08 | 1 | 0 | 0 | 1000 |
| Planted six-leaf hierarchy | Gaussian noise, sigma=.8 | 8 | 0.203 | 23 | 977 |
| Crossed rectangle | Feature | 2 | 0.094 | 12 | 0 |
| Crossed rectangle | Gaussian noise, sigma=.08 | 1 | 0 | 0 | 0 |
| Crossed rectangle | Gaussian noise, sigma=.8 | 3 | 0.023 | 2 | 0 |
| Square with equidistant sides | Feature | 7 | 1.606 | 109 | |
| Square with equidistant sides | Gaussian noise, sigma=.08 | 2 | 1.000 | 0 | |
| Square with equidistant sides | Gaussian noise, sigma=.8 | 14 | 2.487 | 280 | |

The six-leaf target is planted truth. The crossed rectangle's diagonal
pairing is a declared comparator with no truth status. The unperturbed
planted and rectangle configurations each have one reachable Ward
topology. The square has exactly two, the two side pairings. The two
reachable trees differ by rooted RF 4. The floating-point enumerator
agreed at zero and default tie tolerance on these fixtures and reproduced
the two square side pairings independently specified by the Ward costs.

In the square, Gaussian noise at sigma=.08 gave side-pairing frequencies
505/1,000 and 495/1,000 and stayed in the exact family. At sigma=.8, the
same split frequencies were 507/1,000 and 474/1,000, and the diagonal edge
split occurred 19/1,000 times. Those split counts include rooted trees
with different internal clades. Feature resampling gave side-split counts
544/1,000 and 402/1,000 and the diagonal split 54/1,000 times. For the
planted hierarchy, the top edge split survived all feature and Gaussian
perturbations, even when full rooted recovery failed. Its feature
resampling rooted RF histogram was 934, 41 and 25 at RF 0, 2 and 4.

Restoring canonical leaf order before fitting gave the unperturbed tree
in all 3,000 row permutations, and labeled distances were identical.
Fitting the square in shuffled order before restoring output labels
produced 512 copies of the nominal tree and 488 of the other tied tree.
That change reflects order-dependent tie-breaking at identical labeled
distances. The planted hierarchy and rectangle were unchanged either way.

The exact tied family counts compatible Ward outcomes at fixed geometry.
Bootstrap frequencies instead measure stability under a named
perturbation law. Feature resampling changes the Gram matrix by weighting
columns and has an iid-feature interpretation only when the feature units
support it. The six fixture columns include repeated and deterministically
related coordinates, so their resampling frequencies are descriptive.
The measurement-noise bootstrap used fixed, independent Gaussian
perturbations and two prespecified scales. Its higher scale reached
topologies outside the original family. None of these frequencies is a
uniform geometry count or a statistical p-value.

The eight focused geometry tests passed. The analysis checks every
row-permutation distance, the expected tie-family sizes and the planted
unperturbed topology. All 12,000 attempts are represented in the
frequency tables.
