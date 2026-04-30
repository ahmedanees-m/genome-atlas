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
       (ESM-2 dimension, typically 640 for esm2_t30_150M_UR50D) to ``hidden_channels``.  Using a
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
                         embedding size (640 for esm2_t30_150M_UR50D).
        dropout:         Dropout probability applied between layers.
    """

    def __init__(
        self,
        metadata: tuple,
        hidden_channels: int = 128,
        out_channels: int = 128,
        model_type: str = "sage",
        num_layers: int = 2,
        in_features: int = 640,
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
                    # 4-head GAT, heads averaged (concat=False) so output stays
                    # out_ch wide.  Self-loops only for homogeneous edge types
                    # (src == dst); bipartite edges cannot have self-loops.
                    conv_dict[(src, rel, dst)] = GATConv(
                        in_ch, out_ch, heads=4, dropout=dropout,
                        add_self_loops=(src == dst), concat=False,
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
            # Keep a reference to the pre-conv features for two reasons:
            # 1. Source-only node types (e.g. System, which never appears as
            #    an edge destination) are dropped by HeteroConv — pass them
            #    through unchanged to avoid NoneType errors in the next layer.
            # 2. GAT residual: GATConv (bipartite, no self-loops) returns
            #    zeros for any destination node that has *no* incoming edges.
            #    On this graph, 9,991/10,000 Protein nodes have zero System
            #    neighbours, so their layer-1 output is a zero vector and
            #    the ESM-2 diversity loaded by lin_in is completely lost →
            #    embedding collapse (4 unique/10 k).  Adding back the input
            #    features (residual) restores that diversity.  The fix is
            #    applied only for GAT; SAGEConv already concatenates self-
            #    features internally and does not need it.
            prev_x = x_dict
            new_x  = conv(x_dict, edge_index_dict)
            x_dict = {}
            for nt, x in prev_x.items():
                if nt not in new_x:
                    x_dict[nt] = x                         # source-only passthrough
                elif self.model_type == "gat":
                    x_dict[nt] = new_x[nt] + x             # GAT residual: preserve input
                else:
                    x_dict[nt] = new_x[nt]                 # SAGE: no residual needed

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
