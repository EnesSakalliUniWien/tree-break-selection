---
title: Selected Neighborhood Measurability Law
type: analysis
status: draft
updated: 2026-08-10
sources:
  - wiki/analyses/traversal-neighborhood-method-comparison.md
  - wiki/analyses/selected-neighborhood-bottleneck-law.md
  - wiki/sources/sibling-null-prior-interpolation-audit-20260604.md
  - wiki/sources/old-vs-current-method-stack-comparison-20260615.md
  - wiki/sources/selected-neighborhood-spectral-flow-diagnostic-20260616.md
  - wiki/sources/selected-neighborhood-topology-frontier-diagnostic-20260616.md
  - wiki/sources/root-selected-region-overlap-case-family-20260616.md
  - wiki/sources/root-selected-tie-cell-burden-20260616.md
  - wiki/sources/root-selected-mixed-region-law-20260616.md
  - wiki/sources/root-tie-rank-calibration-feasibility-20260616.md
  - wiki/sources/root-tie-rank-selected-null-simulation-pilot-20260616.md
  - wiki/sources/root-tie-rank-null-proposal-frontier-20260616.md
  - wiki/sources/root-tie-rank-generated-neighborhood-replay-20260616.md
  - wiki/sources/root-tie-rank-selected-spectral-generator-targets-20260616.md
  - wiki/sources/root-tie-rank-spectral-lift-parameter-sweep-20260616.md
  - wiki/sources/root-tie-rank-conditioned-coherent-topology-join-20260617.md
  - wiki/sources/root-selected-spectral-tail-law-with-legacy-overlay-20260617.md
  - wiki/sources/root-selected-importance-tail-support-20260617.md
  - wiki/sources/root-tie-rank-target-conditioned-importance-frontier-20260617.md
  - wiki/sources/legacy-internal-spectral-comparison-panel-20260616.md
  - benchmarks/diagnostics/calibration/sibling/nulls/sibling_null_prior_interpolation_audit.py
  - benchmarks/diagnostics/calibration/selected/neighborhood/selected_neighborhood_spectral_flow.py
tags:
  - analysis
  - traversal
  - measurability
  - neighborhood
  - bandwidth
---

# Selected Neighborhood Measurability Law

## Summary

The new traversal logic should be a selected-neighborhood measurability law,
not a fragmentation-penalized score. At each candidate sibling decision, the
method should first ask whether the raw sibling test is measurable under the
current selected-family contract. If not, it may form a diagnostic interpolated
null prior only from valid local ancestor and stable-neighborhood evidence,
with signal-neighborhood attenuation and explicit support checks. If neither
direct measurement nor supported interpolation is available, traversal remains
fail-closed with a bottleneck label.

## Details

Let \(G=(V,E)\) be the selected hierarchy, and let \(v\in V\) be an internal
candidate with children \(c_1,c_2\). Let \(d(x,y)\) denote tree distance, let
\(k_x\) be the local projection dimension or Wald degrees of freedom, and let
\(p_x\in[0,1]\) denote a measured child-edge or sibling-null p-value when such
a value exists.

For a child \(c\), define local evidence sets:

\[
A(c)=\{\text{stopped ancestor edges available for }c\},
\]

\[
T(c)=\{\text{stable tested edges in the local selected neighborhood}\},
\]

\[
S(c)=\{\text{significant signal edges in the local selected neighborhood}\}.
\]

The old bandwidth coordinates are reused as local scales:

\[
\tau_b>0,\qquad \tau_t>0,\qquad \tau_s>0,\qquad h_k>0.
\]

The ancestor and stable-neighborhood weights are:

\[
w_A(c,a)
=
\exp\{-d(c,a)/\tau_b\},
\qquad a\in A(c),
\tag{ancestor support}
\]

\[
w_T(c,t)
=
\exp\{-d(c,t)/\tau_t\}
\exp\left\{
-\frac12
\left(
\frac{\log k_t-\log k_c}{h_k}
\right)^2
\right\},
\qquad t\in T(c).
\tag{stable neighborhood}
\]

Define the support denominator:

\[
D(c)=\sum_{a\in A(c)}w_A(c,a)+\sum_{t\in T(c)}w_T(c,t).
\]

The interpolated local null evidence is defined only when \(D(c)>0\):

\[
\widetilde p_0(c)
=
\frac{
\sum_{a\in A(c)}w_A(c,a)p_a
+\sum_{t\in T(c)}w_T(c,t)p_t
}{
D(c)
}.
\label{eq:interpolated-null-evidence}
\]

Signal-neighborhood evidence does not become calibration support. It only
attenuates the local null prior:

\[
\rho_S(c)
=
\max_{s\in S(c)}
(1-p_s)\exp\{-d(c,s)/\tau_s\},
\qquad
0\le \rho_S(c)\le 1.
\label{eq:signal-attenuation}
\]

The diagnostic child null prior is:

\[
\pi_0(c)=\widetilde p_0(c)(1-\rho_S(c)).
\label{eq:child-null-prior}
\]

There is no mathematical clip in Equation \(\ref{eq:child-null-prior}\). Under
the assumptions \(p_a,p_t,p_s\in[0,1]\), nonnegative weights, positive
bandwidths, and \(D(c)>0\), Equation \(\ref{eq:interpolated-null-evidence}\)
is a convex average and Equation \(\ref{eq:signal-attenuation}\) lies in
\([0,1]\). Therefore \(\pi_0(c)\in[0,1]\). If an implementation obtains a
value outside \([0,1]\) beyond numerical tolerance, the row should be marked
invalid or unsupported rather than made valid by clipping.

For the sibling candidate \(v=(c_1,c_2)\), the old signal-permissive bottleneck
used the smaller child null prior:

\[
\pi_0(v)=\min\{\pi_0(c_1),\pi_0(c_2)\}.
\label{eq:pair-prior}
\]

In the new law, Equation \(\ref{eq:pair-prior}\) is not a production p-value.
It is a diagnostic prior-like quantity that can support a rescue only when the
local support contract also passes.

Define the direct measurability indicator:

\[
M_{\mathrm{dir}}(v)=
\mathbf 1\{
\text{sibling p-value at }v\text{ is directly measurable under the selected}
\text{-family contract}
\}.
\]

Define the interpolation support indicator:

\[
M_{\mathrm{int}}(v)=
\mathbf 1\{
D(c_1)>0,\ D(c_2)>0,\ \tau_b,\tau_t,\tau_s,h_k>0,
\text{ and no required support row is selected-nonnull}
\}.
\]

Define a directed topology coherence predicate \(C(v)\), built from incoming
branch balance, outgoing balance, outgoing edge-norm balance, and the
income/outcome balance product:

\[
C(v)=
\mathbf 1\{
B^{\mathrm{in}}_v\ge b_0,
B^{\mathrm{out}}_v\ge t_0,
E^{\mathrm{out}}_v\ge e_0,
B^{\mathrm{in}}_vB^{\mathrm{out}}_v\ge q_0
\}.
\]

When explicit topology fields are absent, the diagnostic can compute a
selected-tree structural fallback from parent links and descendant leaf counts.
Let \(n_v\) be the descendant leaf count of \(v\), \(p(v)\) its parent, and
\(\operatorname{sib}(v)\) the incoming sibling mass under the same parent. Let
\(\ell(v),r(v)\) be the two largest child descendant masses of \(v\) when
available. Then

\[
\widetilde B^{\mathrm{in}}_v
=
\frac{\min\{n_v,n_{\operatorname{sib}(v)}\}}{n_{p(v)}},
\qquad
\widetilde B^{\mathrm{out}}_v
=
\frac{\min\{n_{\ell(v)},n_{r(v)}\}}{n_v},
\]

and the fallback product is

\[
\widetilde P_v
=
\widetilde B^{\mathrm{in}}_v
\widetilde B^{\mathrm{out}}_v.
\]

This fallback is still only a coherence/localization feature. Root candidates
have no incoming branch under this definition and are labeled as requiring a
separate root-selected topology law; rows whose \(\widetilde P_v\) is below
the balance floor remain fail-closed.

The spectral refinement should be attached to the same coherence predicate,
not used as a direct p-value. Let \(q_x\) be the raw
Marchenko--Pastur-certified signal count at node \(x\), and let
\(U_x^{(r)}\) contain the first \(r\) row-basis eigenvectors for that node.
For a neighboring selected-tree edge \(e=(u,w)\), define

\[
r_e=\min\{q_u,q_w\}.
\]

If \(r_e=0\), the edge has no MP-certified shared mode and is labeled
floor-only. If \(r_e>0\), compare the matched eigenspaces by the singular
values \(\sigma_j(e)\) of

\[
U_u^{(r_e)}U_w^{(r_e)\top}.
\]

This removes arbitrary eigenvector sign flips while still detecting real
rotations of the selected subspace. The normalized chordal distance is

\[
D_U(e)=
\sqrt{
\frac{
r_e-\sum_{j=1}^{r_e}\sigma_j(e)^2
}{r_e}
}.
\]

Eigenvalue flow is measured on log scale:

\[
D_\lambda(e)=
\sqrt{
\frac1{r_e}
\sum_{j=1}^{r_e}
\left(\log\lambda_{u,j}-\log\lambda_{w,j}\right)^2
}.
\]

The diagnostic spectral barrier is

\[
B_{\mathrm{spec}}(e)
=
D_U(e)+D_\lambda(e)+
\frac{|q_u-q_w|}{\max(q_u,q_w,1)}
-\gamma R(e),
\]

where \(R(e)\) is a bounded reward for clear local eigengaps. The spectral
flow affinity is

\[
A_{\mathrm{spec}}(e)=\exp\{-B_{\mathrm{spec}}(e)\}.
\]

Plainly, \(A_{\mathrm{spec}}\) asks whether the same physical mode appears to
flow from leaves toward the parent without rotating away or changing scale too
quickly. The existing bandwidth \(\tau\) variables then act like diffusion
lengths: small \(\tau\) trusts only very local evidence, while large \(\tau\)
allows information to spread farther through the selected tree. The spectral
diagnostic says whether that spread is coherent in the MP-supported modes.

The multiplicity-aware version replaces single-index mode matching with a
block spectral signature. For node \(x\), let

\[
\mathfrak S_x=
\left\{
\left(P_{x,a},m_{x,a},\chi_{x,a},\bar \ell_{x,a}\right)
\right\}_{a=1}^{K_x},
\]

where \(P_{x,a}\) is the projector onto MP block \(a\), \(m_{x,a}\) is its
multiplicity, \(\chi_{x,a}\) is its scale-normalized characteristic
polynomial, and \(\bar\ell_{x,a}\) is its mean log eigenvalue. Across an edge
\(e=(u,w)\), the block transport cost is

\[
D_{\mathrm{mode}}(u,w)=
\min_{\pi}
\frac{1}{\max(K_u,K_w,1)}
\left[
\sum_{(a,b)\in\pi}
\left(
\alpha D_P(P_{u,a},P_{w,b})
+\beta |\bar\ell_{u,a}-\bar\ell_{w,b}|
+\eta D_m(m_{u,a},m_{w,b})
+\zeta D_\chi(\chi_{u,a},\chi_{w,b})
\right)
+\rho N_{\mathrm{unmatched}}
\right].
\]

The corresponding conductance is

\[
c_{uw}
=
\exp\{-D_{\mathrm{mode}}(u,w)\}.
\]

This is the implemented analogue of a connection-Laplacian residual: evidence
is allowed to flow only when MP mode blocks can be transported with low
projector, eigenvalue, multiplicity, and polynomial distortion.

The traversal law is:

\[
\boxed{
\operatorname{action}(v)=
\begin{cases}
\texttt{split},
& M_{\mathrm{dir}}(v)=1,\ p_{\mathrm{sib}}(v)\le \alpha,
\\[3pt]
\texttt{diagnostic\_rescue},
& M_{\mathrm{dir}}(v)=0,\ M_{\mathrm{int}}(v)=1,\ C(v)=1,
\pi_0(v)\le \alpha_{\mathrm{int}},
\\[3pt]
\texttt{fail\_closed},
& \text{otherwise.}
\end{cases}}
\label{eq:selected-neighborhood-action}
\]

The first implementation step is narrower than
Equation \(\ref{eq:selected-neighborhood-action}\): spectral transport is
introduced only as a pass-through support guard. If a closed sibling gate would
otherwise pass through to a descendant split, traversal now requires a
supported measured MP-mode path when `require_mp_blocks=True`. Floor-only,
missing, or unmatched MP evidence is therefore an unmeasured bottleneck rather
than support. The optional non-required diagnostic mode can keep such edges
neutral, but the promoted targeted guard follows the measurability law: no
measured support means no spectral support for pass-through. This realizes the
bottleneck logic without promoting interpolated evidence into a split or
rescue rule.

Fragmentation is not an argument of
Equation \(\ref{eq:selected-neighborhood-action}\). Cluster count, singleton
fraction, and effective cluster count are run-level audit outcomes. They can
show that a law failed, but they should not appear as local inference terms.

The bottleneck label is determined before returning `fail_closed`:

\[
\operatorname{bottleneck}(v)=
\begin{cases}
\texttt{direct\_measurable\_not\_significant},
& M_{\mathrm{dir}}(v)=1,\ p_{\mathrm{sib}}(v)>\alpha,
\\
\texttt{support\_bottleneck},
& M_{\mathrm{dir}}(v)=0,\ D(c_1)D(c_2)=0,
\\
\texttt{bandwidth\_bottleneck},
& \min(\tau_b,\tau_t,\tau_s,h_k)\le0,
\\
\texttt{selection\_bottleneck},
& \text{available neighborhood evidence is selected-nonnull only},
\\
\texttt{topology\_bottleneck},
& M_{\mathrm{int}}(v)=1,\ C(v)=0,
\\
\texttt{interpolated\_null\_not\_strong},
& M_{\mathrm{int}}(v)=1,\ C(v)=1,\ \pi_0(v)>\alpha_{\mathrm{int}}.
\end{cases}
\]

Plainly, the method should split directly when the selected-family contract
gives a valid sibling p-value. It should rescue only when direct measurement is
unavailable, interpolation support is valid, the candidate is coherent in
directed topology, and the interpolated null prior is sufficiently small. It
should otherwise stop and report the mathematical reason.

## Evidence

- [[traversal-neighborhood-method-comparison]] records the old bandwidth
  equations and the distinction between old sibling-null-prior interpolation
  and current strict empirical-null inflation.
- [[selected-neighborhood-bottleneck-law]] records the correction that
  fragmentation is an audit outcome, not a production penalty.
- [[sibling-null-prior-interpolation-audit-20260604]] records why interpolated
  priors are diagnostic unless support excludes selected non-null borrowing.
- [[selected-neighborhood-spectral-flow-diagnostic-20260616]] records that
  MP eigenvalue/eigenvector flow has weak but directionally useful
  signal-vs-selected-null separation on three overlap cases, so it is a
  bottleneck localizer and possible stratum variable rather than a direct
  rescue threshold. Its multiplicity-aware extension adds MP block projectors,
  multiplicities, normalized characteristic polynomials, and mode-transport
  residuals; the current overlap run has almost only singleton blocks, so
  multiplicity and polynomial terms are implemented but not active separators.
- [[selected-neighborhood-topology-frontier-diagnostic-20260616]] records the
  root/non-root threshold version of the same problem. At `tau_s = 20`,
  direct-positive bandwidth reopen counts are `579` selected-null versus `388`
  signal rows per method profile. Root outgoing-balance floors and lowered
  non-root balance-product floors also admit selected-null rows more readily
  than signal rows, so root structural balance and non-root topology frontiers
  remain bottleneck/localization variables until a selected-root margin law and
  stricter spectral/topology support are available.
- [[root-selected-region-overlap-case-family-20260616]] narrows the root-law
  blocker: all seven overlap case-family roots have
  `discrete_tie_cell_geometry_required`, with numerical-zero root-child
  construction margins and hundreds of tied minimum construction merges. The
  missing root law is therefore a discrete selected-region/tie-cell
  conditioning law, not a smooth outgoing-balance threshold.
- [[root-selected-tie-cell-burden-20260616]] quantifies that discrete blocker
  as \(\sum_t \log m_t\) over tied root-child construction choices. The burden
  is large in all seven overlap cases, but it is not monotone with the root
  sibling selected ratio. The selected pair's tie-rank fraction is more
  aligned with the root selected ratio, arguing for a conditional tie-cell law
  using deterministic tie-breaking plus edge/spectral root statistics, not a
  tie-count penalty.
- [[root-selected-mixed-region-law-20260616]] turns that statement into the
  explicit event
  \[
  \mathcal E_{\mathrm{root}}
  =
  \mathcal E_{\mathrm{margin}}
  \cap
  \mathcal E_{\mathrm{tie}}
  \cap
  \mathcal E_{\mathrm{rank}}.
  \]
  All seven overlap roots are `discrete_tie_rank_region` rows with
  `blocked_until_discrete_tie_rank_null_calibrated`; six of seven also show
  bandwidth reopening without root-law support, and no row has hybrid strict
  support. The correct next object is therefore a calibrated conditional law
  for the discrete tie-rank root event, not a bandwidth or topology rescue.
- [[root-tie-rank-calibration-feasibility-20260616]] makes the simulation
  requirement explicit. The seven overlap roots occupy seven conditioning
  strata over tie-rank band, edge-margin band, spectral-ratio band, and
  bandwidth-reopen status. None of the current rows are admissible selected-null
  calibration support. For target \(\alpha=0.01\), plus-one empirical p-value
  resolution requires `99` selected-null roots per stratum, and a `0.25`
  relative tail standard-error target requires `1584` selected-null roots per
  stratum. On the current seven strata this is `693` root simulations for
  resolution only, or `11088` for the stated tail precision.
- [[root-tie-rank-selected-null-simulation-pilot-20260616]] starts that
  simulation path. A one-replicate iid Bernoulli selected-null pilot over the
  seven overlap scales produces `7/7` successful root rows and `0` failures,
  but the generated null rows occupy only `3` null-generated strata and none of
  the seven observed target strata. Thus the immediate blocker is sharper than
  raw null count: the selected-null generator must be conditioned or enriched
  enough to hit the observed high edge-margin, high spectral-ratio, and
  bandwidth-reopen strata before a root tail estimate is meaningful.
- [[root-tie-rank-null-proposal-frontier-20260616]] tests that enrichment
  direction explicitly. In a two-case smoke over `overlap_mod_4c_small` and
  `overlap_mod_6c_med`, iid and column-beta proposal rows still stay at small
  selected ratios, while a two-block tilted proposal produces very large root
  selected ratios and edge margins but low spectral-ratio bands. A sparse
  block-spike proposal reaches `spectral_ratio_gt_4` for one generated root,
  but with low edge margin and low selected ratio. A coupled dense-plus-sparse
  proposal still behaves like the dense action rows: it creates huge selected
  ratios and high edge margins, but the spectral-ratio band remains low. All
  generated roots still miss the observed target strata, and unjoined
  topology-frontier bandwidth is now labeled `bandwidth_reopen_missing` rather
  than measured no-reopen. Thus the missing analytic rule is not merely a
  high-action, high-edge-margin, high-spectral-ratio, or additive
  action-plus-spike generator; it must jointly condition selected tie rank,
  edge margin, spectral ratio, and measured bandwidth-reopen behavior.
- [[root-tie-rank-generated-neighborhood-replay-20260616]] performs that
  generated-matrix replay. All `10` generated proposal matrices replay
  successfully and every generated proposal root now has measured
  topology-frontier evidence. Rebuilding the proposal frontier with those rows
  changes generated root bandwidth from missing to measured: `3/10` generated
  roots reopen at reference bandwidth and `7/10` are measured no-reopen. The
  updated coupling panel has measured-neighborhood evidence for every proposal
  family, but measured-neighborhood coupling still reaches only easier targets
  and not the hard overlap roots. The remaining blocker is therefore no longer
  generated bandwidth coverage; it is the selected spectral-action/tie-rank
  law and external support for the discrete root stratum.
- [[root-tie-rank-selected-spectral-generator-targets-20260616]] turns that
  residual into proposal-family generator requirements. Only the coupled
  edge-spectral and two-block tilt diagnostic families cover `7/7` observed
  targets in the measured high action-edge/tie stratum, and neither reaches
  spectral excess for any target. The coupled family needs median spectral
  lift `2.493244` and maximum `4.296772`; the two-block family needs median
  lift `2.595615` and maximum `4.473194`. The selected-null iid family covers
  `0/7` targets, so the missing object is not a stronger threshold on the
  current selected-null generator. It is a selected-null generator or external
  law that simultaneously occupies the high action-edge/tie measured stratum
  and has the required spectral-excess tail.
- [[root-tie-rank-spectral-lift-parameter-sweep-20260616]] tests whether the
  generator knobs can create that spectral tail before generated-neighborhood
  replay. On `overlap_mod_6c_med`, moderate and stronger coupled
  dense-plus-spike settings cover `7/7` observed targets in root
  action-edge/tie metrics, but both reach `0/7` target spectral excesses. The
  stronger coupled setting increases root selected ratio and edge margin while
  lowering selected spectral-excess log from `0.338471` to `0.244586`. A new
  coherent rank-one spike proposal improves the root screen: the best compact
  setting reaches `2/7` targets and lowers median required spectral lift from
  `2.718312` to `1.739915`. This supports the MP-spike path, but the remaining
  `5/7` misses mean coherent population mode construction is not sufficient by
  itself. A matched-target conditioned coherent spike map, using only selected
  tie-rank/action-edge geometry and no target spectral excess, reaches the same
  `2/7` targets with median residual lift `1.805514` and maximum `2.952730`.
  Thus target conditioning is a useful no-borrowing localization audit, but the
  remaining hard roots still require measured-neighborhood topology replay or a
  sharper selected-root spectral tail.
- [[root-tie-rank-conditioned-coherent-topology-join-20260617]] performs that
  measured-neighborhood topology replay for the conditioned coherent spike
  rows and joins the result back into the spectral generator target panel. All
  seven conditioned coherent rows obtain measured root-frontier evidence, but
  every row is `bandwidth_no_root_reopen`. With matched-target filtering, the
  family still reaches only `2/7` targets and has median residual spectral lift
  `1.805514`. Thus the remaining hard-root blocker is no longer an unrun
  topology replay for this family; it is the selected-root spectral tail and
  external selected-null support.
- [[root-selected-spectral-tail-law-with-legacy-overlay-20260617]] makes that
  final root-tail object explicit. It conditions on the root selected event,
  selected tie-rank \(T\), selected-ratio action \(A\), edge action \(E\),
  measured bandwidth/topology \(B\), and population-law status \(H_u\), while
  leaving \(S_{\mathrm{root}}\) as the tail variable. All seven observed roots
  have `selected_null_support_count = 0` and therefore fail closed. The legacy
  overlay does not contradict this: the full old method leaks one selected-null
  false split on `overlap_mod_4c_small`, while the internal-barycenter spectral
  mode changes MP counts without providing selected-root tail support.
- [[root-selected-importance-tail-support-20260617]] adds the next inference
  path for that blocker. Instead of relabeling tilted proposal rows as ordinary
  null rows, it records \(\log(dP_0/dQ)\) for tilted binary proposals and lets
  the root-tail panel compute an ESS-conservative weighted tail p-value inside
  the same \(T,A,E,B,H_u\) stratum. After generated-neighborhood topology
  replay and joining measured bandwidth evidence, the accumulated importance
  smoke gives same-stratum weighted support for `2/7` observed roots:
  `overlap_extreme_4c` and `overlap_unbal_6c_med` have non-exceeding weighted
  support and conservative p-values `0.5`. The `overlap_unbal_6c_med` stratum
  has six support rows but ESS remains `1.0`, showing weight degeneracy.
  Additional scalar two-block tilts miss the remaining low/mid, mid/high, and
  mid/mid action-edge strata, so the next generator should target
  \(T,A,E,B\) directly rather than continue blind global tilt sweeps.
- [[root-tie-rank-target-conditioned-importance-frontier-20260617]] implements
  the first version of that target-conditioned search, preserving
  likelihood-ratio external-null weights and attaching
  `conditioning_target_case_id` to generated rows. The initial narrow pure
  two-block and coupled smokes on `overlap_heavy_4c_small_feat` and
  `overlap_mod_4c_small` both report `pre_topology_supported_target_count =
  0`. Thus the remaining obstacle is not only target scoping; the current
  Bernoulli tilt families still cannot continuously control the mixed
  action-edge bands needed by the unsupported roots. Adding an unbalanced
  two-block likelihood-ratio proposal gets closer for `overlap_mod_4c_small`,
  but boundary probes still alternate between low/low and high/high
  action-edge bands with no mid/high hit. An accepted-stratum rejection smoke
  retains only rows matching the target pre-topology \(T,A,E\) key, but with
  the current unbalanced proposal it retains `0/20` attempted candidates. The
  next two-factor proposal adds a partially correlated residual root factor to
  test whether parent-space spectral energy can raise edge action without
  forcing selected sibling action high. It also reports `0/8` pre-topology
  hits on `overlap_mod_4c_small`. The missing law therefore needs either a
  sharper selected-root analytic tail or a higher-dimensional proposal that
  directly conditions the selected tie-cell/tie-rank event, not only block
  imbalance, residual factor energy, or post-hoc rejection.
- [[legacy-internal-spectral-comparison-panel-20260616]] records the copied
  old internal-barycenter spectral path as a standard-dispatch diagnostic. It
  greatly increases MP threshold rows and raw MP signal counts, but the
  completed one-replicate overlap partitions are identical to the current
  leaf-only spectral path. This argues against treating internal rows as a
  direct rescue rule for high fragmentation.

## Links

- [[selected-neighborhood-bottleneck-law]]
- [[traversal-neighborhood-method-comparison]]
- [[sibling-null-prior-interpolation-audit-20260604]]
- [[selected-neighborhood-spectral-flow-diagnostic-20260616]]
- [[selected-neighborhood-topology-frontier-diagnostic-20260616]]
- [[root-selected-region-overlap-case-family-20260616]]
- [[root-selected-tie-cell-burden-20260616]]
- [[root-selected-mixed-region-law-20260616]]
- [[root-tie-rank-calibration-feasibility-20260616]]
- [[root-tie-rank-selected-null-simulation-pilot-20260616]]
- [[root-tie-rank-null-proposal-frontier-20260616]]
- [[root-tie-rank-generated-neighborhood-replay-20260616]]
- [[root-tie-rank-selected-spectral-generator-targets-20260616]]
- [[root-tie-rank-spectral-lift-parameter-sweep-20260616]]
- [[root-tie-rank-conditioned-coherent-topology-join-20260617]]
- [[root-selected-spectral-tail-law-with-legacy-overlay-20260617]]
- [[root-selected-importance-tail-support-20260617]]
- [[root-tie-rank-target-conditioned-importance-frontier-20260617]]
- [[legacy-internal-spectral-comparison-panel-20260616]]

## Open Questions

- What threshold \(\alpha_{\mathrm{int}}\) is admissible for an interpolated
  diagnostic prior that is not a calibrated p-value?
- How much MP-supported spectral-flow evidence is required before a local
  bandwidth stratum is identifiable rather than floor-only?
- Which local support count is sufficient for \(M_{\mathrm{int}}(v)=1\) when
  signal-neighborhood evidence is present but excluded from null support?
- Should \(C(v)\) be a hard predicate or a monotone frontier over incoming and
  outgoing topology variables?
