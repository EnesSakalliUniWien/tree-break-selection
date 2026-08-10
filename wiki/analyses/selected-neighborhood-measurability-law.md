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
  - wiki/sources/root-selected-region-overlap-case-family-20260616.md
  - wiki/sources/root-selected-tie-cell-burden-20260616.md
  - wiki/sources/root-selected-mixed-region-law-20260616.md
  - wiki/sources/root-tie-rank-calibration-feasibility-20260616.md
  - wiki/sources/root-tie-rank-selected-null-simulation-pilot-20260616.md
  - wiki/sources/root-tie-rank-null-proposal-frontier-20260616.md
  - wiki/sources/legacy-internal-spectral-comparison-panel-20260616.md
  - benchmarks/diagnostics/calibration/sibling/nulls/sibling_null_prior_interpolation_audit.py
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
- [[root-selected-region-overlap-case-family-20260616]] defines the root-margin
  replay that distinguishes tie-heavy discrete selected regions from smooth
  first-order margin regions.
- [[root-selected-tie-cell-burden-20260616]] represents the discrete component
  as \(\sum_t \log m_t\) over tied construction choices while keeping the
  selected pair's deterministic tie rank separate. Both are conditioning
  coordinates, not calibrated penalties.
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
  The retained contract classifies discrete roots as
  `blocked_until_discrete_tie_rank_null_calibrated`. The correct next object is
  therefore a calibrated conditional law for the discrete tie-rank root event.
- [[root-tie-rank-calibration-feasibility-20260616]] makes the simulation
  requirement explicit over component, tie-rank, edge-margin, and
  spectral-ratio bands. Only explicitly admissible calibration roles count as
  null support.
- [[root-tie-rank-selected-null-simulation-pilot-20260616]] generates selected
  null roots through the canonical root extractor and recomputes the reduced
  feasibility key before target support is summarized.
- [[root-tie-rank-null-proposal-frontier-20260616]] separates proposal-family
  reachability from calibration admissibility and preserves likelihood-ratio
  metadata for supported external-null families.
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
- [[root-selected-region-overlap-case-family-20260616]]
- [[root-selected-tie-cell-burden-20260616]]
- [[root-selected-mixed-region-law-20260616]]
- [[root-tie-rank-calibration-feasibility-20260616]]
- [[root-tie-rank-selected-null-simulation-pilot-20260616]]
- [[root-tie-rank-null-proposal-frontier-20260616]]
- [[legacy-internal-spectral-comparison-panel-20260616]]

## Open Questions

- What threshold \(\alpha_{\mathrm{int}}\) is admissible for an interpolated
  diagnostic prior that is not a calibrated p-value?
- Which local support count is sufficient for \(M_{\mathrm{int}}(v)=1\) when
  signal-neighborhood evidence is present but excluded from null support?
- Should \(C(v)\) be a hard predicate or a monotone frontier over incoming and
  outgoing topology variables?
