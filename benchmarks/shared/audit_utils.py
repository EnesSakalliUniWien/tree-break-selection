"""Utilities for exporting diagnostic audit logs from hierarchical clustering runs.

This module handles the extraction, cleaning, and CSV export of statistical metadata
from tree decompositions, ensuring reports are human-readable and auditable.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import List

import networkx as nx
import numpy as np
import pandas as pd

from benchmarks.shared.result_records import ComputedResultRecord

logger = logging.getLogger(__name__)


def export_decomposition_audit(
    computed_results: List[ComputedResultRecord],
    output_root: Path,
    verbose: bool = False,
) -> None:
    """Export detailed statistical and merge audit logs to CSV files.

    Parameters
    ----------
    computed_results
        List of typed result records from benchmark_cluster_algorithm.
    output_root
        Base directory for benchmark results (e.g. benchmarks/results).
        Logs will be saved in [output_root]/audit/.
    verbose
        If True, print status messages for each saved file.
    """
    audit_dir = output_root / "audit"
    audit_dir.mkdir(parents=True, exist_ok=True)

    for res in computed_results:
        case_id = res.test_case_num
        method = res.method
        tree = res.tree

        # 1. Export Primary Node Stats
        annotations_df = res.annotations
        if annotations_df is not None:
            # ...existing code...
            cols_to_drop = ["distribution"]
            export_df = annotations_df.drop(
                columns=[c for c in cols_to_drop if c in annotations_df.columns]
            ).copy()

            if tree is not None:
                node_labels = nx.get_node_attributes(tree, "label")
                node_parents = {c: p for p, c in tree.edges()}

                # Extract branch lengths from edges
                branch_lengths = {}
                for _p, c, data in tree.edges(data=True):
                    branch_lengths[c] = data.get("branch_length")

                # Enrich with Label, Parent context, and Branch Length
                export_df.insert(
                    0,
                    "node_label",
                    [node_labels.get(idx, idx) for idx in export_df.index],
                )
                export_df["parent_node"] = export_df.index.map(node_parents)
                export_df["parent_label"] = export_df["parent_node"].apply(
                    lambda x: node_labels.get(x, x) if pd.notna(x) else ""
                )
                export_df["branch_length"] = export_df.index.map(branch_lengths)

            # Polishing: Round numeric values for human scanning, but PRESERVE p-values
            p_val_cols = [
                c
                for c in export_df.columns
                if "p_value" in c.lower() or "_bh" in c.lower() or "corrected" in c.lower()
            ]
            numeric_cols = export_df.select_dtypes(include=[np.number]).columns
            round_cols = [c for c in numeric_cols if c not in p_val_cols]

            export_df[round_cols] = export_df[round_cols].round(6)

            stats_file = audit_dir / f"case_{case_id}_{method}_stats.csv"
            export_df.to_csv(stats_file, index=True)
            if verbose:
                print(f"  Node stats saved: {stats_file.name}")
