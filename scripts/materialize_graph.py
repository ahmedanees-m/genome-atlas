"""Materialize DuckDB atlas as NetworkX graph (atlas.gpickle)."""
import argparse
import pickle
from pathlib import Path

import duckdb
import networkx as nx


def materialize(duckdb_path: Path, output_path: Path):
    print("=" * 60)
    print("GENOME-ATLAS Graph Materialization (Step 12)")
    print("=" * 60)

    con = duckdb.connect(str(duckdb_path), read_only=True)

    # Create a directed multi-graph (allows multiple edge types between same nodes)
    G = nx.MultiDiGraph()

    # --- Nodes ---
    node_tables = [
        ("System",    "nodes_system",    ["name", "type", "subtype", "mechanism_bucket"]),
        ("Protein",   "nodes_protein",   ["accession", "sequence", "length", "organism_id", "reviewed", "protein_name"]),
        ("Domain",    "nodes_domain",    ["accession", "name", "source", "mechanism_bucket"]),
        ("Structure", "nodes_structure", ["accession", "source", "method", "resolution_A", "mean_plddt"]),
        ("Mechanism", "nodes_mechanism", ["name", "bucket", "chemistry", "requires_host_repair"]),
        ("Organism",  "nodes_organism",  ["ncbi_taxon_id", "scientific_name", "lineage"]),
        # RNA nodes (added in v0.6.0 — table may not exist in older DBs)
        ("RNA",       "nodes_rna",       ["name", "rna_type", "length_nt"]),
    ]

    total_nodes = 0
    for node_type, table, attrs in node_tables:
        # nodes_rna may not exist in databases built before v0.6.0
        try:
            df = con.execute(f"SELECT * FROM {table}").fetchdf()
        except duckdb.CatalogException:
            print(f"  Nodes [{node_type}]: table '{table}' not found — skipping")
            continue
        for _, row in df.iterrows():
            node_id = f"{node_type}_{int(row['id'])}"
            data = {"node_type": node_type}
            for attr in attrs:
                if attr in row:
                    data[attr] = row[attr]
            G.add_node(node_id, **data)
        total_nodes += len(df)
        print(f"  Nodes [{node_type}]: {len(df):,}")

    # --- Edges ---
    edges_df = con.execute("SELECT * FROM edges").fetchdf()
    for _, row in edges_df.iterrows():
        src = f"{row['source_type']}_{int(row['source_id'])}"
        tgt = f"{row['target_type']}_{int(row['target_id'])}"
        G.add_edge(
            src, tgt,
            edge_type=row["edge_type"],
            weight=float(row["weight"]),
            provenance=row["provenance"],
            confidence=float(row["confidence"]),
        )

    print(f"  Edges: {len(edges_df):,}")
    con.close()

    # --- Stats ---
    print("")
    print("=" * 60)
    print("GRAPH STATISTICS")
    print("=" * 60)
    print(f"  Nodes: {G.number_of_nodes():,}")
    print(f"  Edges: {G.number_of_edges():,}")
    print(f"  Density: {nx.density(G):.6f}")

    if G.number_of_nodes() > 0:
        in_degrees = [d for n, d in G.in_degree()]
        out_degrees = [d for n, d in G.out_degree()]
        print(f"  Avg in-degree: {sum(in_degrees) / len(in_degrees):.2f}")
        print(f"  Avg out-degree: {sum(out_degrees) / len(out_degrees):.2f}")

    # Node type distribution
    print("")
    print("  Node type distribution:")
    type_counts = {}
    for n, d in G.nodes(data=True):
        nt = d.get("node_type", "UNKNOWN")
        type_counts[nt] = type_counts.get(nt, 0) + 1
    for nt, cnt in sorted(type_counts.items()):
        print(f"    {nt}: {cnt:,}")

    # Edge type distribution
    print("")
    print("  Edge type distribution:")
    edge_type_counts = {}
    for u, v, d in G.edges(data=True):
        et = d.get("edge_type", "UNKNOWN")
        edge_type_counts[et] = edge_type_counts.get(et, 0) + 1
    for et, cnt in sorted(edge_type_counts.items()):
        print(f"    {et}: {cnt:,}")

    # --- Save ---
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'wb') as f:
        pickle.dump(G, f)
    sz = output_path.stat().st_size
    print(f"")
    print(f"  Saved: {output_path}")
    print(f"  Size: {sz / (1024 * 1024):.1f} MB")
    print("")
    print("Done!")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--duckdb", type=Path, required=True)
    p.add_argument("--output-gpickle", type=Path, required=True)
    args = p.parse_args()
    materialize(args.duckdb, args.output_gpickle)
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
