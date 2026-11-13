from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING

import networkx as nx

from tree_break_selection.tree.distributions import populate_distributions
from tree_break_selection.tree.elements import Partition, PartitionSet
from tree_break_selection.tree.topology import (
    bottom_up_nodes,
    compute_descendant_leaf_sets,
    get_leaf_values,
    is_leaf,
    lowest_common_ancestor_for_set,
    tree_root,
)

if TYPE_CHECKING:
    import pandas as pd

    from tree_break_selection.tree.feature_space import FeatureSpace


# ============================================================
# 1) PosetTree (NetworkX.DiGraph subclass)
# ============================================================


class PosetTree(nx.DiGraph):
    """Directed tree wrapper that exposes hierarchy operations.

    The class augments ``networkx.DiGraph`` with helpers that make hierarchical
    clustering workflows easier to manage:

    * the root (in-degree 0) is tracked and can be retrieved via :meth:`root`.
    * leaves carry a ``label`` attribute so downstream consumers can recover the
      original sample identifiers.
    * construction modules create a consistent directed representation where
      edges always point from parent to child.
    * utility accessors (:meth:`get_leaves`, :meth:`compute_descendant_sets`)
      provide common tree queries needed by statistical routines.

    Leaves have ``out_degree == 0`` and are expected to carry ``is_leaf=True``.
    """

    def __init__(self, *args, **kwargs):
        """Initialize PosetTree with annotations_df property."""
        super().__init__(*args, **kwargs)
        self.annotations_df: pd.DataFrame | None = None

    # ---------------- Poset helpers ----------------

    def root(self) -> object:
        """Return the cached root node, discovering it if necessary."""
        root = tree_root(self)
        self.graph["root"] = root
        return root

    def get_leaves(self, node: object | None = None, sort: bool = True) -> list[object]:
        """Return leaf labels, globally or under ``node``; the method form of
        ``topology.get_leaf_values``."""
        return get_leaf_values(self, node=node, sort=sort)

    def compute_descendant_sets(self, use_labels: bool = True) -> dict[object, frozenset]:
        """Map each node to the set of leaf labels under it; the method form of
        ``topology.compute_descendant_leaf_sets``.

        Parameters
        ----------
        use_labels
            When ``True`` (default), map to stored ``label`` values; otherwise use
            internal node identifiers.

        Returns
        -------
        dict[str, frozenset]
            Dictionary whose keys are node ids and whose values are the descendant leaf
            labels/ids as a frozenset.
        """
        return compute_descendant_leaf_sets(self, use_labels=use_labels)

    # ---------------- Set based identity ----------------

    def initialize_split_indices(self, encoding: dict[object, int] | None = None) -> None:
        """Give every node its leaf set as ``split_indices`` and index the nodes by it.

        A node's ``split_indices`` is a ``Partition`` of its descendant leaves, so the same
        leaf set has the same identity in every tree that shares ``encoding`` (leaf label to
        bit position), whatever its node id is. By default leaves are numbered in insertion
        order (``L0``, ``L1``, ... for linkage trees). Call again after changing the topology.
        """
        leaves = [node for node in self.nodes if is_leaf(self, node)]
        if encoding is None:
            encoding = {self.nodes[leaf]["label"]: index for index, leaf in enumerate(leaves)}
        if len({encoding[self.nodes[leaf]["label"]] for leaf in leaves}) != len(leaves):
            raise ValueError("Leaf labels must map to distinct encoding positions.")

        split_index: dict[Partition, object] = {}
        for node in bottom_up_nodes(self):
            if is_leaf(self, node):
                bitmask = 1 << encoding[self.nodes[node]["label"]]
            else:
                bitmask = 0
                for child in self.successors(node):
                    bitmask |= self.nodes[child]["split_indices"].bitmask
            split = Partition.from_bitmask(bitmask, encoding)
            self.nodes[node]["split_indices"] = split
            split_index[split] = node
        self.graph["taxa_encoding"] = encoding
        self.graph["split_index"] = split_index

    def to_splits(self, with_leaves: bool = False) -> PartitionSet[Partition]:
        """Return the leaf sets of the internal nodes (and of the leaves on request)."""
        if "split_index" not in self.graph:
            self.initialize_split_indices()
        return PartitionSet(
            {
                split
                for split, node in self.graph["split_index"].items()
                if with_leaves or not is_leaf(self, node)
            },
            encoding=self.graph["taxa_encoding"],
        )

    def find_node_by_split(self, target_split: Partition | tuple[int, ...]) -> object | None:
        """Return the node whose leaf set is ``target_split``, or ``None`` if it does not exist."""
        if "split_index" not in self.graph:
            self.initialize_split_indices()
        encoding = self.graph["taxa_encoding"]
        if not isinstance(target_split, Partition):
            target_split = Partition(tuple(target_split), encoding)
        elif target_split.encoding is not encoding and target_split.encoding != encoding:
            raise ValueError("Cannot search for a split with a different encoding.")
        return self.graph["split_index"].get(target_split)

    def find_lca_for_set(self, nodes: Iterable[object]) -> object:
        """Find the lowest common ancestor for a collection of nodes.

        Iteratively applies the two-node LCA function to find the LCA for the set.

        Parameters
        ----------
        nodes
            An iterable of node identifiers.

        Returns
        -------
        str
            The node id of the lowest common ancestor for the set.
        """
        return lowest_common_ancestor_for_set(self, nodes)

    def populate_node_distributions(
        self,
        leaf_data: "pd.DataFrame",
        *,
        feature_space: "FeatureSpace | None" = None,
    ) -> None:
        """Populate tree nodes with distributions and build stats DataFrame.

        Populates each node with:
        - distribution: empirical subtree barycenter, computed as the
          leaf-count-weighted mean of leaf/child distributions
        - leaf_count: number of descendant leaves

        Node distributions are stored as flat raw-coordinate vectors. The
        optional feature-space contract supplies Bernoulli/categorical/continuous block
        structure for downstream covariance code.

        Parameters
        ----------
        leaf_data
            DataFrame where rows are leaf labels and columns are raw feature coordinates.

        Notes
        -----
        Results are stored in ``tree.annotations_df`` for later access.
        """
        import pandas as pd

        populate_distributions(
            self,
            leaf_data,
            feature_space=feature_space,
        )
        node_records = []
        for node_id in self.nodes():
            node_attrs = self.nodes[node_id]
            node_records.append(
                {
                    "node_id": node_id,
                    "distribution": node_attrs["distribution"],
                    "leaf_count": node_attrs["leaf_count"],
                    "is_leaf": node_attrs["is_leaf"],
                }
            )
        self.annotations_df = pd.DataFrame.from_records(node_records).set_index(
            "node_id", drop=True
        )
