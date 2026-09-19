# GO analysis plotting library research — 2026-09-19

Requested scope: deduplicate plotting in the adaptive diffusion cosine subspace
GO-IC pipeline, give each plot type a separate Python module, and separate
DataFrame input/output from computation. No dependency installation was requested.

## Documentation retrieval

Context7 library IDs resolved and queried:

- `/websites/matplotlib_stable`: reusable explicit Figure/Axes plotting helpers.
- `/websites/seaborn_pydata`: DataFrame heatmaps, supplied Axes, labels and color limits.
- `/scipy/scipy`: vectorized binary entropy from SciPy special functions.

Official pages supporting the comparison:

- https://matplotlib.org/stable/users/explain/quick_start.html
- https://matplotlib.org/stable/users/explain/figure/api_interfaces.html
- https://seaborn.pydata.org/generated/seaborn.heatmap.html
- https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html
- https://github.com/scipy/scipy/blob/main/doc/source/tutorial/special.rst

The Matplotlib, Seaborn and SciPy dendrogram pages were also opened directly.
These are current documentation pages, not a claim that their current release
versions are installed. Local versions checked: Matplotlib 3.10.7, SciPy 1.16.3,
pandas 2.3.3, NumPy 2.3.4, scikit-learn 1.7.2, pypdf 6.12.1.

## Findings and decisions

| Candidate | What it replaces or supports | Decision |
| --- | --- | --- |
| Matplotlib explicit Axes API | Repeated chart rendering can be one drawing function reused by figure constructors and report export. | Adopt the documented composition pattern using the installed dependency. Figure constructors take arrays/DataFrames and return figures; export owns file saving. |
| SciPy `dendrogram(..., ax=...)` | Tree geometry and leaf order are already supplied by SciPy. Its returned leaf order supports a cluster strip without custom tree traversal. | Keep the existing dependency and centralize the repeated SciPy call in `draw_dendrogram`. |
| Seaborn `heatmap` | Handles DataFrame labels, colorbars and configurable color limits on a supplied Axes. | Viable future replacement for heatmap boilerplate. Do not add it in this refactor: it is not a declared project dependency, and its pcolormesh rendering would change the existing imshow output. Retain term selection as separate analysis logic. |
| SciPy `xlogy` and `xlog1py` | Can express binary entropy with endpoint-safe special functions. | Research candidate only. Entropy is outside plotting deduplication; preserve its current numerical and invalid-input behavior. |

The decision to retain the current rendering backend is an implementation
judgment based on dependency scope and output compatibility. It is not a claim
that Seaborn cannot render these data. No new library is needed for deduplication.
