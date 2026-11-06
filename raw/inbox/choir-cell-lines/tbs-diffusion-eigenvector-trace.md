# Left and right diffusion eigenvectors and the root split

On 2026-11-04 we recovered the 30 retained eigenvectors for the saved graph on
48,879 cells. The four trees share this graph and spectrum. Their embeddings
and linkage rules determine how the modes contribute to the root.

## Left and right eigenvectors

Let \(W\) be the symmetric kernel and \(D\) its diagonal matrix of row sums.
The saved symmetric coordinates are \(U\Lambda^3\), where

\[
S=D^{-1/2}WD^{-1/2},\qquad SU=U\Lambda,\qquad U^\top U=I.
\]

For the random-walk operator \(P=D^{-1}W\), the corresponding right and left
eigenvectors, both stored as columns, are

\[
R=D^{-1/2}U,\qquad L=D^{1/2}U,\qquad
PR=R\Lambda,\qquad P^\top L=L\Lambda,\qquad L^\top R=I.
\]

The Markov trees use \(R\Lambda^3\). The left vectors give the dual weights
for the same eigenvalues. The maximum relative eigen-residuals were
\(1.80\times10^{-14}\) for the right vectors and \(1.82\times10^{-14}\) for
the left vectors. The maximum entrywise error in \(L^\top R-I\) was
\(1.55\times10^{-14}\). [Coifman and Lafon, Appendix A](https://www.math.ucdavis.edu/~strohmer/courses/270/diffusion_maps.pdf)

## Component modes and the global cutoff

The eigenvalue 1 has multiplicity three. For a component \(C\) with degree
volume \(v_C=\sum_{i\in C}d_i\), one component-supported dual pair is

\[
r_C=\frac{\mathbf{1}_C}{\sqrt{v_C}},\qquad
l_C=\frac{D\mathbf{1}_C}{\sqrt{v_C}}.
\]

The right vector is constant on its component, and the left vector weights
cells by degree. The saved block solver orders these modes as A204, IM95 and
main, up to sign. The stationary subspace is invariant to the choice of
basis within this repeated eigenvalue.

| Graph component | Cells | Retained nonstationary modes | Largest nonstationary eigenvalue in the block |
| --- | ---: | ---: | ---: |
| Main | 48,444 | 27 | 0.999971375 |
| A204 | 250 | 0 | 0.936761595 |
| IM95 | 185 | 0 | 0.863919346 |

The global cutoff is 0.999051795. A204 and IM95 each have only a stationary
coordinate among the 30 retained modes. Their maximum within-component range
in any saved Markov coordinate is at most \(5.76\times10^{-16}\). The saved
symmetric coordinates vary within those blocks through the square root of
degree. The omitted local modes contain the remaining within-component
diffusion geometry.

## How the modes enter the root split

For Ward, the contributions are the squared coordinate differences between
the two root-child centroids. For average linkage, they are the mean squared
cross-child coordinate differences. Average linkage itself takes the mean
Euclidean distance. We independently recomputed all four saved root heights
using their actual linkage rules.

| Tree | Smaller root group | Stationary share of squared measure | Largest contribution | Saved root height | Height without stationary modes, same split |
| --- | --- | ---: | --- | ---: | ---: |
| Symmetric / average | ASPC1 singleton | 0.4453% | Mode 26, 66.07% | 0.121647 | 0.121385 |
| Symmetric / Ward | A204 | 99.9997% | Stationary subspace | 1.397068 | 0.002284 |
| Markov / average | OAW28 | 0.5346% | Mode 15, 75.29% | 0.022063 | 0.022009 |
| Markov / Ward | NCIH727 | 0.00443% | Mode 17, 95.02% | 0.376794 | 0.376786 |

The last column evaluates the same two children in the reduced coordinates.
It describes the geometry of each recorded split. Modes 15, 17 and 26 all
live in the main component. Their eigenvalues are 0.999626186, 0.999550944
and 0.999242901, respectively.

## A single hierarchy

A component-aware single tree can allocate nonstationary modes within each
component and constrain its cells to form a clade. The upper merges can then
use a declared dissimilarity in the common PCA space. This gives the hierarchy
one root and retains the local variation currently lost from A204 and IM95.
The spectral trace evaluates the saved trees. The proposed merge rule and its
TBS gate behavior remain to be tested.

## Evidence

- `tbs-diffusion-eigenvector-trace.csv` lists all eigenvalues, component support
  and left and right residuals.
- `tbs-diffusion-root-mode-contributions.csv` lists every mode's contribution
  and the coordinate means of the two children in linkage order.
- `tbs-diffusion-eigenvector-trace-evidence.json` records component volumes,
  residuals, local eigenvalues, root summaries and verified input hashes.
- The local run `results/analyses/choir_diffusion_eigenvector_trace_20261104_104311/`
  contains `markov_eigenvectors.npz`. Its `left` and `right` arrays each have
  48,879 cell rows and 30 mode columns. `cell_id` and `eigenvalues` identify
  their rows and paired columns. The calculation is in
  `tbs-diffusion-eigenvector-trace-script.py.txt`.
