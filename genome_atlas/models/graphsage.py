"""Heterogeneous GNN (GraphSAGE / GAT) and link predictor for GENOME-ATLAS.

Used by train_gnn.py for link-prediction benchmarking on the atlas graph.
Both models follow the same HeteroConv + linear-projection architecture so
that the only difference is the per-edge-type aggregation kernel.
"""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import HeteroConv, SAGEConv, GATConv


class HeteroGNN(nn.Module):
    """Two-layer heterogeneous GNN for the atlas knowledge graph.

    Architecture
    ------------
    1. Linear input projection: maps every node type from ``in_features``
       (ESM-2 dimension, typically 480) to ``hidden_channels``.  Using a
       shared projection keeps the interface clean when non-Protein nodes
       have zero features.
    2. ``num_layers`` HeteroConv layers, each wrapping either SAGEConv
       (``model_type='sage'``) or GATConv (``model_type='gat'``) for every
       edge type declared in the graph metadata.
    3. ReLU activation + dropout between layers; no activation on the final
       layer so downstream predictors can apply their own non-linearity.

    Args:
        metadata:        Tuple (node_types, edge_types) returned by
                         ``HeteroData.metadata()``.
        hidden_channels: Width of the hidden layers.
        out_channels:    Dimensionality of the final node embeddings.
        model_type:      ``'sage'`` for GraphSAGE, ``'gat'`` for GAT.
        num_layers:      Number of message-passing layers (default 2).
        in_features:     Input feature dimension; should match the ESM-2
                         embedding size (default 480 for esm2_t30_150M).
        dropout:         Dropout probability applied between layers.
    """

    def __init__(
        self,
        metadata: tuple,
        hidden_channels: int = 128,
        out_channels: int = 128,
        model_type: str = "sage",
        num_layers: int = 2,
        in_features: int = 480,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        self.model_type = model_type
        self.num_layers = num_layers
        self.dropout = dropout

        node_types, edge_types = metadata

        # Per-node-type input projection (shared weight per node type).
        # Projecting from in_features → hidden_channels avoids lazy-init
        # issues with HeteroConv when some node types have zero features.
        self.lin_in = nn.ModuleDict(
            {nt: nn.Linear(in_features, hidden_channels) for nt in node_types}
        )

        self.convs = nn.ModuleList()
        for layer_idx in range(num_layers):
            in_ch  = hidden_channels
            out_ch = out_channels if layer_idx == num_layers - 1 else hidden_channels

            conv_dict: dict = {}
            for (src, rel, dst) in edge_types:
                if model_type == "gat":
                    # Single-head GAT; no self-loops (bipartite edges present)
                    conv_dict[(src, rel, dst)] = GATConv(
                        in_ch, out_ch, heads=1, dropout=dropout,
                        add_self_loops=False, concat=False,
                    )
                else:  # default: sage
                    conv_dict[(src, rel, dst)] = SAGEConv(in_ch, out_ch)

            self.convs.append(HeteroConv(conv_dict, aggr="mean"))

    def forward(
        self,
        x_dict: dict[str, torch.Tensor],
        edge_index_dict: dict[tuple, torch.Tensor],
    ) -> dict[str, torch.Tensor]:
        # Input projection — only project node types present in the graph
        x_dict = {
            nt: F.relu(self.lin_in[nt](x))
            for nt, x in x_dict.items()
            if nt in self.lin_in
        }

        for i, conv in enumerate(self.convs):
            x_dict = conv(x_dict, edge_index_dict)
            is_last = i == self.num_layers - 1
            if not is_last:
                x_dict = {k: F.relu(v) for k, v in x_dict.items()}
            x_dict = {
                k: F.dropout(v, p=self.dropout, training=self.training)
                for k, v in x_dict.items()
            }
        return x_dict


class LinkPredictor(nn.Module):
    """MLP link predictor: concatenates two embeddings → P(edge exists).

    Architecture: Linear(2·d → d) → ReLU → Linear(d → 1) → Sigmoid.
    Returns a 1-D tensor of probabilities, one per input pair.

    Args:
        in_channels: Dimensionality of each node embedding (out_channels of
                     HeteroGNN).
    """

    def __init__(self, in_channels: int) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(in_channels * 2, in_channels),
            nn.ReLU(),
            nn.Linear(in_channels, 1),
            nn.Sigmoid(),
        )

    def forward(self, z_src: torch.Tensor, z_dst: torch.Tensor) -> torch.Tensor:
        """
        Args:
            z_src: (N, d) source node embeddings.
            z_dst: (N, d) destination node embeddings.

        Returns:
            (N,) tensor of edge existence probabilities.
        """
        z = torch.cat([z_src, z_dst], dim=-1)
        return self.net(z).squeeze(-1)
