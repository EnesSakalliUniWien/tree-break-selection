"""Sweep edge/sibling alphas for selected KAK/cosine diagnostic lenses.

The main adaptive-diffusion clustering is kept fixed as the context tree.
This diagnostic changes only the lens tree alpha decisions, then measures how
each lens refines or fragments the main contexts.
"""

from __future__ import annotations

import argparse
import math
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import product
from pathlib import Path

from tree_break_selection.plot.backend import configure_matplotlib_backend

configure_matplotlib_backend()

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
)

from benchmarks.diagnostics.spectral.adaptive_cosine.adaptive_cosine_kak_benchmark_probe import (
    SCHEMA_VERSION,
    sibling_method_counts,
)
from benchmarks.diagnostics.spectral.adaptive_cosine.adaptive_cosine_kak_matrix_probe import (
    load_matrix,
)
from tree_break_selection.hierarchy_analysis.cluster_assignments import (
    build_sample_cluster_assignments,
)
from tree_break_selection.hierarchy_analysis.decomposition.gates.orchestrator import (
    run_gate_annotation_pipeline,
)
from tree_break_selection.hierarchy_analysis.statistics.sibling_divergence.inflation_correction.types.inflation_model import (
    DEFAULT_INTERNAL_SUPPORT_THRESHOLDS,
    CalibrationSupportThresholds,
)
from tree_break_selection.hierarchy_analysis.tree_decomposition import TreeDecomposition
from tree_break_selection.space_separation import (
    SpectralBlock,
    adaptive_spectral_blocks,
    block_adaptive_diffusion_geometry,
    block_diffusion_geometry,
    coordinates_for_block,
    cosine_eigendecomposition,
    weight_feature_matrix,
)
from tree_break_selection.tree.construction import build_tree
from tree_break_selection.tree.poset_tree import PosetTree

SCHEMA = f"{SCHEMA_VERSION}/kak_lens_alpha_sweep"


@dataclass(frozen=True)
class LensSpec:
    family: str
    weighting: str
    block_name: str

    @property
    def lens_id(self) -> str:
        return f"{self.family}__{self.weighting}__{self.block_name}"


DEFAULT_LENSES = (
    LensSpec("raw_kak", "tfidf", "adaptive_modes_14_30"),
    LensSpec("raw_kak", "raw", "adaptive_modes_10_15"),
    LensSpec("sep_fixed_diffusion", "tfidf", "adaptive_modes_14_30"),
    LensSpec("sep_fixed_diffusion", "raw", "adaptive_modes_02_05"),
    LensSpec("sep_adaptive_diffusion", "raw", "adaptive_modes_02_05"),
    LensSpec("sep_adaptive_diffusion", "tfidf", "adaptive_modes_06_09"),
)


def parse_float_list(value: str) -> list[float]:
    values = [float(part.strip()) for part in value.split(",") if part.strip()]
    if not values:
        raise argparse.ArgumentTypeError("expected at least one comma-separated float")
    return values


def parse_linkage_methods(value: str) -> list[str]:
    methods = [part.strip() for part in value.split(",") if part.strip()]
    if not methods:
        raise argparse.ArgumentTypeError("expected at least one linkage method")
    allowed = {"average", "complete", "ward"}
    invalid = sorted(set(methods) - allowed)
    if invalid:
        raise argparse.ArgumentTypeError(f"unsupported linkage method(s): {', '.join(invalid)}")
    return methods


def parse_alpha_pairs(value: str) -> list[tuple[float, float]]:
    pairs: list[tuple[float, float]] = []
    for raw_pair in value.split(","):
        raw_pair = raw_pair.strip()
        if not raw_pair:
            continue
        pieces = raw_pair.split(":")
        if len(pieces) != 2:
            raise argparse.ArgumentTypeError(
                "alpha pairs must be edge:sibling entries, e.g. 0.001:0.01"
            )
        pairs.append((float(pieces[0]), float(pieces[1])))
    if not pairs:
        raise argparse.ArgumentTypeError("expected at least one alpha pair")
    return pairs


def build_alpha_pairs(
    edge_alphas: list[float],
    sibling_alphas: list[float],
) -> list[tuple[float, float]]:
    return [(float(edge), float(sibling)) for edge, sibling in product(edge_alphas, sibling_alphas)]


def parse_lens(value: str) -> LensSpec:
    parts = value.split(":")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError(
            "lens must be family:weighting:block_name, e.g. raw_kak:tfidf:adaptive_modes_14_30"
        )
    family, weighting, block_name = parts
    if family not in {"raw_kak", "sep_fixed_diffusion", "sep_adaptive_diffusion"}:
        raise argparse.ArgumentTypeError(f"unknown lens family {family!r}")
    if weighting not in {"raw", "tfidf"}:
        raise argparse.ArgumentTypeError(f"unknown weighting {weighting!r}")
    return LensSpec(family=family, weighting=weighting, block_name=block_name)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sweep edge/sibling alphas for selected KAK/cosine lenses."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/feature_matrices/feature_matrix_julia_GOCC_GOBP_GOMF_combined.tsv"),
    )
    parser.add_argument(
        "--main-assignments",
        type=Path,
        default=Path(
            "benchmarks/results/00_current_20260427_blob_analysis/20260427_blob_analysis/"
            "09_full_data_method_reference_comparison_20260610/tbs_diffusion_adaptive/"
            "cluster_assignments.csv"
        ),
    )
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument(
        "--edge-alphas",
        type=parse_float_list,
        default=parse_float_list("0.0003,0.001,0.003"),
    )
    parser.add_argument(
        "--sibling-alphas",
        type=parse_float_list,
        default=parse_float_list("0.003,0.01,0.03"),
    )
    parser.add_argument(
        "--alpha-pairs",
        type=parse_alpha_pairs,
        default=None,
        help="Optional comma-separated edge:sibling pairs. Overrides the grid product.",
    )
    parser.add_argument(
        "--tree-linkage-methods",
        type=parse_linkage_methods,
        default=parse_linkage_methods("average"),
        help="Comma-separated tree linkage methods over Euclidean lens coordinates: average,ward.",
    )
    parser.add_argument("--max-rank", type=int, default=80)
    parser.add_argument("--min-segment-length", type=int, default=4)
    parser.add_argument("--max-segments", type=int, default=8)
    parser.add_argument("--diffusion-k-neighbors", type=int, default=15)
    parser.add_argument("--diffusion-time", type=int, default=3)
    parser.add_argument("--diffusion-components", type=int, default=30)
    parser.add_argument("--adaptive-bandwidth-type", default="-1/(d+2)")
    parser.add_argument("--adaptive-epsilon", default="median")
    parser.add_argument("--adaptive-metric", choices=("euclidean",), default="euclidean")
    parser.add_argument(
        "--lens",
        action="append",
        type=parse_lens,
        default=None,
        help="Repeatable family:weighting:block_name override. Defaults to the current top lenses.",
    )
    parser.add_argument(
        "--enforce-internal-support-thresholds",
        action="store_true",
        help="Fail closed when internal empirical-null support is below thresholds.",
    )
    return parser.parse_args()


def default_output_dir(input_path: Path) -> Path:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%SZ")
    return (
        Path("benchmarks/results/diagnostics") / f"kak_lens_alpha_sweep_{input_path.stem}_{stamp}"
    )


def load_main_labels(path: Path, index: pd.Index) -> np.ndarray:
    labels = pd.read_csv(path, index_col=0)
    if "cluster_id" not in labels.columns:
        raise ValueError(f"{path} must contain a cluster_id column")
    return labels.reindex(index)["cluster_id"].astype(int).to_numpy()


def find_block(blocks: list[SpectralBlock], name: str) -> SpectralBlock:
    for block in blocks:
        if block.block_name == name:
            return block
    available = ", ".join(block.block_name for block in blocks)
    raise ValueError(f"Block {name!r} not found. Available: {available}")


def build_lens_coordinates(
    *,
    eigvals: np.ndarray,
    eigvecs: np.ndarray,
    block: SpectralBlock,
    spec: LensSpec,
    args: argparse.Namespace,
) -> tuple[np.ndarray, dict[str, object]]:
    coords = coordinates_for_block(eigvals, eigvecs, block)
    if spec.family == "raw_kak":
        return coords, {"tree_metric": "raw_kak_block_euclidean"}
    elif spec.family == "sep_fixed_diffusion":
        geometry = block_diffusion_geometry(
            coords,
            k_neighbors=args.diffusion_k_neighbors,
            diffusion_time=args.diffusion_time,
            n_components=args.diffusion_components,
        )
        return geometry.coordinates, geometry.metadata
    elif spec.family == "sep_adaptive_diffusion":
        geometry = block_adaptive_diffusion_geometry(
            coords,
            k_neighbors=args.diffusion_k_neighbors,
            diffusion_time=args.diffusion_time,
            n_components=args.diffusion_components,
            metric=args.adaptive_metric,
            bandwidth_type=args.adaptive_bandwidth_type,
            epsilon=args.adaptive_epsilon,
        )
        return geometry.coordinates, geometry.metadata

    raise ValueError(f"unknown lens family {spec.family!r}")


def build_lens_tree(
    *,
    data: pd.DataFrame,
    lens_coordinates: np.ndarray,
    spec: LensSpec,
    tree_linkage_method: str,
    metadata: dict[str, object],
) -> tuple[PosetTree, dict[str, object]]:
    coords = np.asarray(lens_coordinates, dtype=float)
    if coords.ndim != 2 or coords.shape[0] < 3 or not np.isfinite(coords).all():
        raise ValueError(f"invalid coordinates for {spec.lens_id}")

    if tree_linkage_method in {"average", "complete"}:
        tree_metric = "euclidean_condensed"
    elif tree_linkage_method == "ward":
        tree_metric = "euclidean_observations"
    else:
        raise ValueError(f"unsupported tree_linkage_method={tree_linkage_method!r}")
    distances = pdist(coords, metric="euclidean")
    if not np.isfinite(distances).all() or np.allclose(distances, 0.0):
        raise ValueError(f"degenerate distances for {spec.lens_id}")
    # Ward on Euclidean observations equals Ward on their condensed Euclidean distances.
    tree = build_tree(data, distances, linkage_method=tree_linkage_method).tree
    tree.populate_node_distributions(data)
    return tree, {
        **metadata,
        "tree_linkage_method": tree_linkage_method,
        "tree_metric": tree_metric,
    }


def labels_from_assignments(assignments: pd.DataFrame, index: pd.Index) -> np.ndarray:
    return assignments.reindex(index)["cluster_id"].astype(int).to_numpy()


def weighted_purity(labels: np.ndarray, reference: np.ndarray) -> float:
    table = pd.crosstab(pd.Series(labels, name="labels"), pd.Series(reference, name="reference"))
    if table.empty:
        return math.nan
    return float(table.max(axis=1).sum() / table.to_numpy().sum())


def split_counts(labels: np.ndarray, reference: np.ndarray) -> dict[str, int]:
    frame = pd.DataFrame({"lens": labels, "main": reference})
    cluster_sizes = frame["lens"].value_counts()
    cross_context = 0
    cross_context_medium = 0
    for lens_id, group in frame.groupby("lens", sort=False):
        size = int(cluster_sizes.loc[lens_id])
        if size < 3:
            continue
        purity = float(group["main"].value_counts().max() / size)
        if purity < 0.70:
            cross_context += 1
            if size >= 10:
                cross_context_medium += 1

    main_split_ge2 = 0
    main_split_ge3 = 0
    max_subcontexts = 0
    for _, group in frame.groupby("main", sort=False):
        sub_sizes = group["lens"].value_counts()
        n_sub = int((sub_sizes >= 3).sum())
        max_subcontexts = max(max_subcontexts, n_sub)
        if n_sub >= 2:
            main_split_ge2 += 1
        if n_sub >= 3:
            main_split_ge3 += 1

    return {
        "n_cross_context_clusters_ge3_purity_lt70": int(cross_context),
        "n_cross_context_medium_clusters_purity_lt70": int(cross_context_medium),
        "n_main_contexts_split_ge2_subcontexts": int(main_split_ge2),
        "n_main_contexts_strongly_split_ge3_subcontexts": int(main_split_ge3),
        "max_subcontexts_in_one_main_context": int(max_subcontexts),
    }


def score_assignments(
    assignments: pd.DataFrame, main_labels: np.ndarray, index: pd.Index
) -> dict[str, object]:
    labels = labels_from_assignments(assignments, index)
    sizes = pd.Series(labels).value_counts()
    split = split_counts(labels, main_labels)
    lens_purity = weighted_purity(labels, main_labels)
    main_purity = weighted_purity(main_labels, labels)
    singleton_genes = int((sizes == 1).sum())
    score = (
        lens_purity
        + 0.25 * min(split["n_main_contexts_strongly_split_ge3_subcontexts"], 5) / 5.0
        + 0.25 * min(split["max_subcontexts_in_one_main_context"], 30) / 30.0
        - singleton_genes / len(labels)
        - 0.002 * split["n_cross_context_clusters_ge3_purity_lt70"]
    )
    return {
        "n_clusters": int(len(sizes)),
        "singletons": singleton_genes,
        "gene_singleton_fraction": float(singleton_genes / len(labels)),
        "median_cluster_size": float(sizes.median()),
        "largest_cluster_size": int(sizes.max()),
        "largest_cluster_fraction": float(sizes.max() / len(labels)),
        "lens_cluster_purity_to_main_context": float(lens_purity),
        "main_context_purity_to_lens": float(main_purity),
        "nmi_to_main_context": float(normalized_mutual_info_score(main_labels, labels)),
        "ari_to_main_context": float(adjusted_rand_score(main_labels, labels)),
        "alpha_lens_score": float(score),
        **split,
    }


def run_alpha_grid(
    *,
    data: pd.DataFrame,
    main_labels: np.ndarray,
    tree: PosetTree,
    alpha_pairs: list[tuple[float, float]],
    enforce_internal_support_thresholds: bool,
    internal_support_thresholds: CalibrationSupportThresholds = DEFAULT_INTERNAL_SUPPORT_THRESHOLDS,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for edge_alpha, sibling_alpha in alpha_pairs:
        start = time.perf_counter()
        try:
            gate_bundle = run_gate_annotation_pipeline(
                tree,
                tree.annotations_df.copy(),
                edge_alpha=edge_alpha,
                sibling_alpha=sibling_alpha,
                leaf_data=data,
                enforce_internal_support_thresholds=enforce_internal_support_thresholds,
                internal_support_thresholds=internal_support_thresholds,
            )
            decomposition = TreeDecomposition(
                tree=tree,
                gate_annotation_bundle=gate_bundle,
            ).decompose_tree()
            assignments = build_sample_cluster_assignments(decomposition).loc[data.index]
            rows.append(
                {
                    "edge_alpha": float(edge_alpha),
                    "sibling_alpha": float(sibling_alpha),
                    "status": "ok",
                    "runtime_sec": float(time.perf_counter() - start),
                    "sibling_test_method_counts": sibling_method_counts(tree.annotations_df),
                    **score_assignments(assignments, main_labels, data.index),
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "edge_alpha": float(edge_alpha),
                    "sibling_alpha": float(sibling_alpha),
                    "status": "failed_gate",
                    "runtime_sec": float(time.perf_counter() - start),
                    "error": repr(exc),
                }
            )
    return rows


def write_plots(summary: pd.DataFrame, output_dir: Path) -> None:
    import matplotlib.pyplot as plt

    ok = summary[summary["status"].eq("ok")].copy()
    if ok.empty:
        return

    fig, ax = plt.subplots(figsize=(10.8, 7.2))
    for lens_id, group in ok.groupby("lens_tree_id", sort=False):
        ax.plot(
            group["n_clusters"],
            group["lens_cluster_purity_to_main_context"],
            marker="o",
            linestyle="",
            label=lens_id,
            alpha=0.78,
        )
    ax.set_xlabel("Lens clusters")
    ax.set_ylabel("Lens purity to main adaptive context")
    ax.set_title("KAK/cosine lens alpha sweep: fragmentation vs purity")
    ax.grid(alpha=0.25)
    ax.legend(fontsize=7, loc="best")
    fig.tight_layout()
    fig.savefig(output_dir / "alpha_sweep_fragmentation_vs_purity.png", dpi=180)
    plt.close(fig)

    ranked = ok.sort_values("alpha_lens_score", ascending=False).head(20).copy()
    labels = [
        f"{row.lens_tree_id}\ne={row.edge_alpha:g}, s={row.sibling_alpha:g}"
        for row in ranked.itertuples(index=False)
    ]
    fig_height = max(5.0, 0.42 * len(ranked) + 1.8)
    fig, ax = plt.subplots(figsize=(12.0, fig_height))
    ax.barh(np.arange(len(ranked)), ranked["alpha_lens_score"], color="#4c78a8", alpha=0.86)
    ax.set_yticks(np.arange(len(ranked)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("Diagnostic lens score")
    ax.set_title("Top KAK/cosine lens alpha settings")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "alpha_sweep_top_lens_scores.png", dpi=180)
    plt.close(fig)


def write_report(summary: pd.DataFrame, output_dir: Path) -> None:
    ok = summary[summary["status"].eq("ok")].copy()
    lines = [
        "# KAK/Cosine Lens Alpha Sweep",
        "",
        "Main context tree is fixed to full adaptive diffusion. This sweep changes only",
        "the selected KAK/cosine lens tree edge and sibling alphas.",
        "",
        "## Status Counts",
        "",
        summary["status"].value_counts(dropna=False).to_string(),
        "",
    ]
    if not ok.empty:
        top_cols = [
            "lens_id",
            "tree_linkage_method",
            "edge_alpha",
            "sibling_alpha",
            "n_clusters",
            "gene_singleton_fraction",
            "lens_cluster_purity_to_main_context",
            "main_context_purity_to_lens",
            "n_main_contexts_strongly_split_ge3_subcontexts",
            "max_subcontexts_in_one_main_context",
            "n_cross_context_clusters_ge3_purity_lt70",
            "alpha_lens_score",
        ]
        top = ok.sort_values("alpha_lens_score", ascending=False).head(20)
        lines.extend(
            [
                "## Top Alpha Settings",
                "",
                top[top_cols].to_markdown(index=False),
                "",
                "## Best Setting Per Lens",
                "",
                ok.sort_values("alpha_lens_score", ascending=False)
                .groupby("lens_id", as_index=False)
                .head(1)[top_cols]
                .to_markdown(index=False),
                "",
            ]
        )
    (output_dir / "kak_lens_alpha_sweep_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    output_dir = args.output_dir or default_output_dir(args.input)
    output_dir.mkdir(parents=True, exist_ok=True)

    data = load_matrix(args.input)
    main_labels = load_main_labels(args.main_assignments, data.index)
    lens_specs = tuple(args.lens) if args.lens else DEFAULT_LENSES
    alpha_pairs = (
        list(args.alpha_pairs)
        if args.alpha_pairs is not None
        else build_alpha_pairs(args.edge_alphas, args.sibling_alphas)
    )

    eigensystems: dict[
        str, tuple[np.ndarray, np.ndarray, list[SpectralBlock], dict[str, object]]
    ] = {}
    rows: list[dict[str, object]] = []
    for spec in lens_specs:
        print(f"[lens] {spec.lens_id}", flush=True)
        if spec.weighting not in eigensystems:
            values = weight_feature_matrix(data, spec.weighting)
            eigvals, eigvecs = cosine_eigendecomposition(values, args.max_rank)
            blocks, diagnostics = adaptive_spectral_blocks(
                eigvals,
                min_segment_length=args.min_segment_length,
                max_segments=args.max_segments,
            )
            eigensystems[spec.weighting] = (eigvals, eigvecs, blocks, diagnostics)
        eigvals, eigvecs, blocks, diagnostics = eigensystems[spec.weighting]
        block = find_block(blocks, spec.block_name)
        try:
            lens_coordinates, coordinate_metadata = build_lens_coordinates(
                eigvals=eigvals,
                eigvecs=eigvecs,
                block=block,
                spec=spec,
                args=args,
            )
            for tree_linkage_method in args.tree_linkage_methods:
                print(f"[lens] {spec.lens_id} / {tree_linkage_method}", flush=True)
                tree, tree_metadata = build_lens_tree(
                    data=data,
                    lens_coordinates=lens_coordinates,
                    spec=spec,
                    tree_linkage_method=tree_linkage_method,
                    metadata=coordinate_metadata,
                )
                alpha_rows = run_alpha_grid(
                    data=data,
                    main_labels=main_labels,
                    tree=tree,
                    alpha_pairs=alpha_pairs,
                    enforce_internal_support_thresholds=args.enforce_internal_support_thresholds,
                )
                for row in alpha_rows:
                    rows.append(
                        {
                            "schema_version": SCHEMA,
                            "input_path": str(args.input),
                            "main_assignments": str(args.main_assignments),
                            "n_samples": int(data.shape[0]),
                            "n_features": int(data.shape[1]),
                            "lens_family": spec.family,
                            "weighting": spec.weighting,
                            "block_name": spec.block_name,
                            "lens_id": spec.lens_id,
                            "tree_linkage_method": tree_linkage_method,
                            "lens_tree_id": f"{spec.lens_id}__{tree_linkage_method}",
                            "block_start": int(block.block_start),
                            "block_end": int(block.block_end),
                            "subspace_dimensions": int(block.block_end - block.block_start + 1),
                            "segmentation_bic": diagnostics.get("bic", math.nan),
                            "tree_metadata": repr(tree_metadata),
                            **row,
                        }
                    )
        except Exception as exc:
            rows.append(
                {
                    "schema_version": SCHEMA,
                    "input_path": str(args.input),
                    "main_assignments": str(args.main_assignments),
                    "n_samples": int(data.shape[0]),
                    "n_features": int(data.shape[1]),
                    "lens_family": spec.family,
                    "weighting": spec.weighting,
                    "block_name": spec.block_name,
                    "lens_id": spec.lens_id,
                    "tree_linkage_method": ",".join(args.tree_linkage_methods),
                    "lens_tree_id": spec.lens_id,
                    "status": "failed_tree",
                    "error": repr(exc),
                }
            )
        pd.DataFrame(rows).to_csv(output_dir / "kak_lens_alpha_sweep_summary.csv", index=False)

    summary = pd.DataFrame(rows)
    summary.to_csv(output_dir / "kak_lens_alpha_sweep_summary.csv", index=False)
    if not summary.empty and "alpha_lens_score" in summary.columns:
        summary[summary["status"].eq("ok")].sort_values("alpha_lens_score", ascending=False).head(
            30
        ).to_csv(output_dir / "kak_lens_alpha_sweep_top_rows.csv", index=False)
    write_plots(summary, output_dir)
    write_report(summary, output_dir)
    print(f"Wrote KAK/cosine lens alpha sweep: {output_dir}", flush=True)


if __name__ == "__main__":
    main()
