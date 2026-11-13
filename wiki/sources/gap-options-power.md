---
title: Power of Four Gap Options on the Refused Layouts
type: source
status: draft
updated: 2026-11-11
sources:
  - raw/inbox/choir-cell-lines/gap-options-power-20261111.md
  - raw/inbox/choir-cell-lines/gap-options-power-scripts-20261111.py.txt
tags:
  - calibration
  - sibling-gate
  - simulation
---

# Power of Four Gap Options on the Refused Layouts

## Summary

A one-direction simulation compared four additions to `rho` (two-Gaussian
mixture, dip statistic, valley depth, shape of the two sides) on the layouts
that the Gaussian best-flat-cut reference refuses. No option accepts the
connected `cancer_pc` layout. Dip and valley accept the mirrored distant group
and keep the level on skewed and heavy-tailed populations; the mixture and the
shape reject every skewed single population. No option would change the
benchmark result.

## Key Points

- Connected layout (256 cells at 0, 64 at -3 and 3): at most 9 of 1,000
  samples accepted at 0.01 by any option. Along its axis it has valleys at
  0.999 of the neighbouring mode, skewness 0 and kurtosis 3.000.
- 64 cells at offset 3 next to 256: `rho` 0, valley 4, dip 3 of 1,000; mixture
  1,000 of 1,000.
- Mirrored group at offset 20: dip, valley and mixture 1,000 of 1,000; `rho`
  and shape 0.
- One gamma population (shape 4): mixture and best-cut shape reject 2,000 of
  2,000; dip and valley 15 of 2,000. Student t with 5 degrees of freedom:
  mixture 1,403, dip 4, valley 0.
- Equal sizes at offset 3: valley 855, dip 664, `rho` 643 of 1,000. 64 cells
  at offset 4: `rho` 932, valley 546, dip 247.
- Benchmark panel: the connected cases stay unsolved; the disconnected cases
  are already correct through the pass-through traversal.
- Limits: one direction, no tree and no complete fit; the valley's band width
  is a chosen value; the dip is an own implementation; a combined rule of
  `rho` or gap was not simulated.

## Evidence

`raw/inbox/choir-cell-lines/gap-options-power-20261111.md` holds the setup, the
three result tables and the limits; protocol, script and summary are in
`raw/inbox/choir-cell-lines/gap-options-power-scripts-20261111.py.txt`.

## Links

- [[best-flat-cut-reference-benchmark]]
- [[gaussian-best-flat-cut-reference]]
- [[single-population-split-bound]]
- [[ceiling-rule-literature-comparison]]
