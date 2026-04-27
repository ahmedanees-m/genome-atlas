"""Graph view selection for GENOME-ATLAS.

Provides ``get_graph()`` - a high-level entry point that wraps
``build_pyg_hetero`` with a ``graph_view`` parameter so callers can choose
between the compact ML-training graph (``'primary'``) and the full annotated
graph (``'full'``).

Primary view (default, backward-compatible with v0.7.0)
--------------------------------------------------------
Contains only the four edge types used for GNN training::

    Protein   --HAS_DOMAIN-->    Domain
    Structure --STRUCTURE_OF-->  Protein
    System    --USES_MECHANISM--> Mechanism
    System    --HAS_PROTEIN-->   Protein

Full view
---------
Includes all primary edges plus three secondary edge types that were present
in v0.6.0 but omitted from the v0.7.0 training graph for ML-benchmark clarity::

    System   --HAS_RNA-->    RNA      (derived from rna_components in foundational_systems.yaml)
    Protein  --PART_OF-->    System   (reverse of HAS_PROTEIN)
    Protein  --SIMILAR_TO--> Protein  (cosine similarity >= threshold on ESM-2 embeddings)

Downstream packages that need the full annotation graph should call
``get_graph(graph_view='full')``.
"""
from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:  # pragma: no cover
    import torch
    from torch_geometric.data import HeteroData

from genome_atlas.graph.build import add_train_val_test_split, build_pyg_hetero

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Edge types that belong to the primary (ML-training) view
_PRIMARY_EDGE_TYPES: frozenset[tuple[str, str, str]] = frozenset(
    {
        ("Protein", "HAS_DOMAIN", "Domain"),
        ("Structure", "STRUCTURE_OF", "Protein"),
        ("System", "USES_MECHANISM", "Mechanism"),
        ("System", "HAS_PROTEIN", "Protein"),
    }
)

GraphView = Literal["primary", "full"]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def get_graph(
    gpickle_path: Path,
    esm_emb_path: Path,
    *,
    graph_view: GraphView = "primary",
    yaml_path: Path | None = None,
    similarity_threshold: float = 0.90,
    add_split: bool = True,
    train_ratio: float = 0.8,
    val_ratio: float = 0.1,
    seed: int = 42,
) -> HeteroData:
    """Load the GENOME-ATLAS as a PyG HeteroData object.

    Parameters
    ----------
    gpickle_path:
        Path to an atlas gpickle (full or training graph).
    esm_emb_path:
        Path to the ESM-2 parquet (columns: ``accession``, ``embedding``).
    graph_view:
        ``'primary'`` (default) - 4-edge ML-training graph, backward-compatible
        with v0.7.0.  ``'full'`` - primary edges plus HAS_RNA, PART_OF, and
        cosine-similarity SIMILAR_TO edges between Protein nodes.
    yaml_path:
        Path to ``foundational_systems.yaml``.  Required for deriving HAS_RNA
        edges in the ``'full'`` view.  Defaults to the bundled data file
        inside the installed ``genome_atlas`` package.
    similarity_threshold:
        Minimum cosine similarity (0-1) for a SIMILAR_TO edge between two
        Protein nodes.  Only used when ``graph_view='full'``.  Default 0.90.
    add_split:
        If ``True`` (default) apply an 80/10/10 train/val/test mask via
        :func:`~genome_atlas.graph.build.add_train_val_test_split`.
    train_ratio, val_ratio, seed:
        Forwarded to :func:`~genome_atlas.graph.build.add_train_val_test_split`
        when ``add_split=True``.

    Returns
    -------
    HeteroData
        PyG heterogeneous graph ready for training or inference.

    Raises
    ------
    ValueError
        If ``graph_view`` is not ``'primary'`` or ``'full'``.

    Examples
    --------
    >>> # ML training - compact, backward-compatible
    >>> data = get_graph(gpickle_path, esm_emb_path, graph_view='primary')

    >>> # Full annotation graph (all edge types)
    >>> data = get_graph(gpickle_path, esm_emb_path, graph_view='full')
    """
    # Build base graph (includes all edge types present in the gpickle)
    data = build_pyg_hetero(gpickle_path, esm_emb_path)

    if graph_view == "primary":
        data = _filter_primary(data)

    elif graph_view == "full":
        # Retain any extra edges already in the gpickle (e.g. from older
        # versions that stored HAS_RNA / PART_OF / SIMILAR_TO as real edges),
        # then derive any that are still missing.
        data = _add_full_view_edges(
            data,
            yaml_path=yaml_path,
            similarity_threshold=similarity_threshold,
        )

    else:
        raise ValueError(
            f"graph_view must be 'primary' or 'full', got {graph_view!r}"
        )

    if add_split:
        data = add_train_val_test_split(
            data, train_ratio=train_ratio, val_ratio=val_ratio, seed=seed
        )

    return data


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _filter_primary(data: HeteroData) -> HeteroData:  # type: ignore[name-defined]
    """Remove edge types that are not in the primary ML-training view."""
    to_remove = [et for et in data.edge_types if tuple(et) not in _PRIMARY_EDGE_TYPES]
    for et in to_remove:
        del data[et]
    return data


def _add_full_view_edges(
    data: HeteroData,  # type: ignore[name-defined]
    yaml_path: Path | None,
    similarity_threshold: float,
) -> HeteroData:  # type: ignore[name-defined]
    """Derive and add secondary edge types to the full-view graph."""
    import torch  # noqa: PLC0415

    # ------------------------------------------------------------------ #
    # PART_OF  (Protein -> System, reverse of HAS_PROTEIN)
    # ------------------------------------------------------------------ #
    has_protein_key = ("System", "HAS_PROTEIN", "Protein")
    part_of_key = ("Protein", "PART_OF", "System")

    if part_of_key not in data.edge_types and has_protein_key in data.edge_types:
        hp_ei = data[has_protein_key].edge_index  # [2, E]: [sys_idx, prot_idx]
        # Reverse: src=protein, dst=system
        part_of_ei = torch.stack([hp_ei[1], hp_ei[0]], dim=0)
        data[part_of_key].edge_index = part_of_ei
        print(f"  Derived PART_OF edges: {part_of_ei.size(1)}")

    # ------------------------------------------------------------------ #
    # HAS_RNA  (System -> RNA, from YAML rna_components)
    # ------------------------------------------------------------------ #
    has_rna_key = ("System", "HAS_RNA", "RNA")

    if has_rna_key not in data.edge_types:
        if yaml_path is None:
            # Fall back to bundled data file in the installed package
            yaml_path = (
                Path(__file__).parent.parent / "data" / "foundational_systems.yaml"
            )

        if yaml_path.exists():
            _add_has_rna_from_yaml(data, yaml_path)
        else:
            print(
                f"  [warn] yaml_path not found; skipping HAS_RNA derivation: {yaml_path}"
            )

    # ------------------------------------------------------------------ #
    # SIMILAR_TO  (Protein -> Protein, cosine similarity)
    # ------------------------------------------------------------------ #
    similar_to_key = ("Protein", "SIMILAR_TO", "Protein")

    if similar_to_key not in data.edge_types and "Protein" in data.node_types:
        prot_feats = data["Protein"].x  # [N, emb_dim]
        # Only compute if we have non-trivial embeddings
        non_zero_rows = (prot_feats.abs().sum(dim=1) > 0).sum().item()
        if non_zero_rows > 1:
            srcs, dsts = _compute_similarity_edges(prot_feats, similarity_threshold)
            if srcs.numel() > 0:
                data[similar_to_key].edge_index = torch.stack([srcs, dsts], dim=0)
                print(
                    f"  Derived SIMILAR_TO edges: {srcs.numel()} "
                    f"(cosine >= {similarity_threshold}, {non_zero_rows} non-zero proteins)"
                )
            else:
                print(
                    f"  SIMILAR_TO: no pairs above threshold {similarity_threshold}"
                )

    return data


def _add_has_rna_from_yaml(data: HeteroData, yaml_path: Path) -> None:  # type: ignore[name-defined]
    """Populate (System, HAS_RNA, RNA) edges from foundational_systems.yaml."""
    import torch  # noqa: PLC0415

    try:
        import yaml  # pyyaml
    except ImportError:
        print("  [warn] PyYAML not installed; skipping HAS_RNA derivation")
        return

    with open(yaml_path, "r", encoding="utf-8") as f:
        doc = yaml.safe_load(f)

    # Build lookup: node_id -> positional index within the node type
    systems_in_graph: dict[str, int] = {}
    if "System" in data.node_types:
        for idx, nid in enumerate(data["System"].node_ids):
            systems_in_graph[nid] = idx

    rnas_in_graph: dict[str, int] = {}
    if "RNA" in data.node_types:
        for idx, nid in enumerate(data["RNA"].node_ids):
            rnas_in_graph[nid] = idx

    if not systems_in_graph or not rnas_in_graph:
        print("  [warn] System or RNA nodes missing from graph; skipping HAS_RNA")
        return

    sys_indices: list[int] = []
    rna_indices: list[int] = []

    for sys_entry in doc.get("systems", []):
        sys_name = sys_entry.get("name", "")
        if sys_name not in systems_in_graph:
            continue
        for rna_name in sys_entry.get("rna_components", []):
            if rna_name in rnas_in_graph:
                sys_indices.append(systems_in_graph[sys_name])
                rna_indices.append(rnas_in_graph[rna_name])

    if sys_indices:
        ei = torch.tensor([sys_indices, rna_indices], dtype=torch.long)
        data[("System", "HAS_RNA", "RNA")].edge_index = ei
        print(f"  Derived HAS_RNA edges: {len(sys_indices)}")
    else:
        print("  HAS_RNA: no matching System->RNA pairs found in YAML")


def _compute_similarity_edges(
    feats: torch.Tensor,  # type: ignore[name-defined]
    threshold: float,
    batch_size: int = 500,
) -> tuple[torch.Tensor, torch.Tensor]:  # type: ignore[name-defined]
    """Return (src, dst) index pairs where cosine similarity >= threshold.

    Uses batched matrix multiplication to avoid OOM on large protein sets.
    Self-loops are excluded.  The result is directed (both (i,j) and (j,i)
    are included for symmetric pairs).

    Args:
        feats:      [N, D] float tensor of protein embeddings.
        threshold:  Minimum cosine similarity to emit an edge.
        batch_size: Row-batch size for the similarity computation.

    Returns:
        Tuple of (src_indices, dst_indices) 1-D long tensors.
    """
    import torch  # noqa: PLC0415

    n = feats.size(0)
    # L2-normalise rows (zero rows stay zero; clamp avoids divide-by-zero)
    norms = feats.norm(dim=1, keepdim=True).clamp(min=1e-8)
    normed = feats / norms

    src_list: list[torch.Tensor] = []
    dst_list: list[torch.Tensor] = []

    for start in range(0, n, batch_size):
        end = min(start + batch_size, n)
        block = normed[start:end]          # [B, D]
        sim = torch.mm(block, normed.t())  # [B, N]

        # Set self-loop entries to -1 (minimum cosine similarity) so they
        # are always excluded regardless of the threshold value.
        for i in range(end - start):
            sim[i, start + i] = -1.0

        mask = sim >= threshold
        rows, cols = mask.nonzero(as_tuple=True)
        src_list.append(rows + start)
        dst_list.append(cols)

    if not src_list:
        return torch.empty(0, dtype=torch.long), torch.empty(0, dtype=torch.long)

    src = torch.cat(src_list)
    dst = torch.cat(dst_list)
    return src, dst
