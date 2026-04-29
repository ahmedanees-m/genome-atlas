"""Convert NetworkX atlas gpickle + ESM-2 embeddings to PyG HeteroData.

Used by train_gnn.py for heterogeneous GNN training.
"""
from __future__ import annotations

import pickle
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from torch_geometric.data import HeteroData


def build_pyg_hetero(gpickle_path: Path, esm_emb_path: Path) -> HeteroData:
    """Build a PyG HeteroData object from the atlas gpickle and ESM-2 embeddings.

    Node features
    -------------
    - **Protein** nodes: ``emb_dim``-dimensional ESM-2 mean-pooled embeddings
      (640-dim for ``esm2_t30_150M_UR50D``).  Rows with no matching accession
      default to zeros.
    - **All other node types** (Domain, Structure, Mechanism, RNA, System):
      zero vectors of the same ``emb_dim``.  This is intentional:

      * *Domain* and *Mechanism* are purely topological — they function as
        destination-only aggregation targets and acquire meaningful
        representations through message passing from Protein nodes.
      * *Structure* nodes are sources in ``STRUCTURE_OF`` and ``SIMILAR_TO``
        edges.  In layer 1 they contribute a constant (bias-only) signal; by
        layer 2 ``SIMILAR_TO`` updates give them structural neighbourhood
        context.  This is acceptable for the primary Protein→Domain benchmark.
      * *RNA* nodes are destination-only (``System→RNA``); their zero
        features never propagate to other node types.  RNA sequence features
        will be incorporated in a future release via a nucleic-acid language
        model (e.g. RNA-FM).
      * *System* nodes connect sparsely (9 ``HAS_PROTEIN`` edges) — their
        contribution to Protein representations is negligible.

    The embedding dimension is inferred automatically from the ESM-2 parquet,
    so the function works with any ESM-2 variant (35 M → 640 M).

    Edge representation: one ``(src_type, edge_label, dst_type)`` entry per
    directed edge type found in the NetworkX graph.

    Each node type stores a ``node_ids`` attribute (list of NetworkX node ID
    strings, in the same order as the feature matrix rows) so that embeddings
    can be matched back to graph nodes after training.

    Args:
        gpickle_path: Path to atlas_train.gpickle (NetworkX DiGraph,
            isolated nodes already removed by create_train_gpickle.py).
        esm_emb_path: Path to ESM-2 parquet with columns
            ``['accession', 'embedding']``.

    Returns:
        PyG HeteroData ready for training (no train/val/test masks yet;
        call :func:`add_train_val_test_split` before training).
    """
    # ------------------------------------------------------------------ #
    # Load graph
    # ------------------------------------------------------------------ #
    print(f"Loading graph: {gpickle_path}")
    with open(str(gpickle_path), "rb") as f:
        G = pickle.load(f)
    print(f"  {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")

    # ------------------------------------------------------------------ #
    # Load ESM-2 embeddings
    # ------------------------------------------------------------------ #
    print(f"Loading ESM-2 embeddings: {esm_emb_path}")
    esm_df = pd.read_parquet(str(esm_emb_path))
    # Embedding column may be a list or a numpy array; normalise to 1-D ndarray
    first_emb = esm_df["embedding"].iloc[0]
    if not isinstance(first_emb, np.ndarray):
        first_emb = np.array(first_emb, dtype=np.float32)
    emb_dim = len(first_emb)
    esm_map: dict[str, np.ndarray] = {}
    for _, row in esm_df.iterrows():
        emb = row["embedding"]
        if not isinstance(emb, np.ndarray):
            emb = np.array(emb, dtype=np.float32)
        esm_map[row["accession"]] = emb
    print(f"  {len(esm_map):,} ESM-2 embeddings, dim={emb_dim}")

    # ------------------------------------------------------------------ #
    # Group nodes by type; build integer index within each type
    # ------------------------------------------------------------------ #
    nodes_by_type: dict[str, list[str]] = defaultdict(list)
    for node_id, attrs in G.nodes(data=True):
        nt = attrs.get("node_type", "Unknown")
        nodes_by_type[nt].append(node_id)

    # Sort for determinism
    for nt in nodes_by_type:
        nodes_by_type[nt].sort()

    node_to_idx: dict[str, int] = {}
    for nt, node_list in nodes_by_type.items():
        for i, nid in enumerate(node_list):
            node_to_idx[nid] = i

    # ------------------------------------------------------------------ #
    # Build feature matrices
    # ------------------------------------------------------------------ #
    data = HeteroData()

    for nt, node_list in nodes_by_type.items():
        n = len(node_list)
        feats = np.zeros((n, emb_dim), dtype=np.float32)
        if nt == "Protein":
            for i, nid in enumerate(node_list):
                acc = G.nodes[nid].get("accession", "")
                if acc in esm_map:
                    feats[i] = esm_map[acc]
        data[nt].x = torch.from_numpy(feats)
        # Store node ID strings so embeddings can be mapped back after training
        data[nt].node_ids = node_list

    print(f"  Node types: {list(nodes_by_type.keys())}")

    # ------------------------------------------------------------------ #
    # Build edge index tensors grouped by (src_type, edge_label, dst_type)
    # ------------------------------------------------------------------ #
    edges_by_type: dict[tuple, tuple[list, list]] = defaultdict(lambda: ([], []))

    for u, v, edge_attrs in G.edges(data=True):
        et  = edge_attrs.get("edge_type", "UNKNOWN")
        src_type = G.nodes[u].get("node_type", "Unknown")
        dst_type = G.nodes[v].get("node_type", "Unknown")
        key = (src_type, et, dst_type)
        if u in node_to_idx and v in node_to_idx:
            edges_by_type[key][0].append(node_to_idx[u])
            edges_by_type[key][1].append(node_to_idx[v])

    for (src_type, et, dst_type), (srcs, dsts) in edges_by_type.items():
        edge_index = torch.tensor([srcs, dsts], dtype=torch.long)
        data[(src_type, et, dst_type)].edge_index = edge_index

    print(f"  Edge types: {list(edges_by_type.keys())}")
    return data


def add_train_val_test_split(
    data: HeteroData,
    train_ratio: float = 0.8,
    val_ratio: float = 0.1,
    seed: int = 42,
) -> HeteroData:
    """Add train / val / test boolean masks to every edge type in HeteroData.

    Splits are applied per edge type independently, using a fixed random seed
    for reproducibility across runs.

    Args:
        data:        PyG HeteroData (from build_pyg_hetero).
        train_ratio: Fraction of edges used for training (default 0.80).
        val_ratio:   Fraction used for validation (default 0.10).
                     Remainder goes to test.
        seed:        NumPy RNG seed for determinism.

    Returns:
        The same HeteroData with .train_mask / .val_mask / .test_mask on
        each edge-type storage object.
    """
    rng = np.random.default_rng(seed)

    for edge_type in data.edge_types:
        n_edges = data[edge_type].edge_index.size(1)

        if n_edges == 0:
            data[edge_type].train_mask = torch.zeros(0, dtype=torch.bool)
            data[edge_type].val_mask   = torch.zeros(0, dtype=torch.bool)
            data[edge_type].test_mask  = torch.zeros(0, dtype=torch.bool)
            continue

        perm    = rng.permutation(n_edges)
        n_train = int(n_edges * train_ratio)
        n_val   = int(n_edges * val_ratio)

        train_mask = torch.zeros(n_edges, dtype=torch.bool)
        val_mask   = torch.zeros(n_edges, dtype=torch.bool)
        test_mask  = torch.zeros(n_edges, dtype=torch.bool)

        train_mask[perm[:n_train]]                 = True
        val_mask  [perm[n_train:n_train + n_val]]  = True
        test_mask [perm[n_train + n_val:]]         = True

        data[edge_type].train_mask = train_mask
        data[edge_type].val_mask   = val_mask
        data[edge_type].test_mask  = test_mask

    return data
