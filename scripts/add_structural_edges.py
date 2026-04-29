"""Add Structure-SIMILAR_TO-Structure edges from Foldseek all-vs-all results.

Reads foldseek/results.m8 (produced by run_foldseek.sh) and inserts
SIMILAR_TO edges between Structure nodes for pairs with TM-score >= 0.5.

PDB IDs in the m8 file include the .pdb.gz suffix; we strip extensions
and upper-case to match the accession column in nodes_structure.

Keeps at most --max-per-node top-scoring partners per structure to avoid
O(N^2) edge explosion. Removes any previous SIMILAR_TO edges before inserting
(safe for re-runs after Foldseek parameter changes).

Usage (on VM):
    python3 scripts/add_structural_edges.py \
        --duckdb   /home/.../atlas.duckdb \
        --m8       /home/.../foldseek/results.m8 \
        --min-tm   0.5 \
        --max-per-node 10
"""
import argparse
from pathlib import Path
from collections import defaultdict

import duckdb
import pandas as pd


def strip_pdb_id(raw: str) -> str:
    """Convert '3zso.pdb.gz' or '3zso.cif.gz' -> '3ZSO'."""
    stem = Path(raw).name
    for ext in (".pdb.gz", ".pdb", ".cif.gz", ".cif", ".ent.gz", ".ent"):
        if stem.lower().endswith(ext):
            stem = stem[: -len(ext)]
            break
    return stem.upper()


def main(duckdb_path, m8_path, min_tm, max_per_node):
    print("=" * 60)
    print("Structural Similarity Edge Addition")
    print("=" * 60)

    con = duckdb.connect(str(duckdb_path))

    # Build accession → node_id lookup for Structure nodes
    struct_lookup = {
        r[0].upper(): r[1]
        for r in con.execute("SELECT accession, id FROM nodes_structure").fetchall()
    }
    print(f"  Structure nodes in DB  : {len(struct_lookup)}")

    # ------------------------------------------------------------------ #
    # Read m8 results
    # ------------------------------------------------------------------ #
    print(f"  Reading {m8_path} ...")
    # Column layout from run_foldseek.sh --format-output (last col = qtmscore)
    cols = ["query", "target", "fident", "alnlen", "mismatch", "gapopen",
            "qstart", "qend", "tstart", "tend", "evalue", "bits", "tmscore"]
    df = pd.read_csv(str(m8_path), sep="\t", header=None, names=cols,
                     dtype={"tmscore": float})
    # Foldseek outputs qtmscore as the last column; rename for clarity
    if "tmscore" not in df.columns and "qtmscore" in df.columns:
        df = df.rename(columns={"qtmscore": "tmscore"})
    print(f"  Raw hits               : {len(df):,}")

    # Normalise IDs
    df["q_id"] = df["query"].apply(strip_pdb_id)
    df["t_id"] = df["target"].apply(strip_pdb_id)

    # Remove self-hits and duplicates (keep q < t to avoid storing both directions)
    df = df[df["q_id"] != df["t_id"]]
    df = df[df["tmscore"] >= min_tm]
    print(f"  After TM >= {min_tm} filter: {len(df):,}")

    # Map to DB node IDs
    df = df[df["q_id"].isin(struct_lookup) & df["t_id"].isin(struct_lookup)].copy()
    df["src_node"] = df["q_id"].map(struct_lookup)
    df["tgt_node"] = df["t_id"].map(struct_lookup)
    print(f"  After ID mapping       : {len(df):,}")

    # Keep top-N per query node (highest TM-score)
    df = (
        df.sort_values("tmscore", ascending=False)
          .groupby("src_node")
          .head(max_per_node)
          .reset_index(drop=True)
    )
    print(f"  After top-{max_per_node} per node    : {len(df):,}")

    # ------------------------------------------------------------------ #
    # Remove any existing SIMILAR_TO edges (idempotent re-run)
    # ------------------------------------------------------------------ #
    old = con.execute(
        "SELECT COUNT(*) FROM edges WHERE edge_type='SIMILAR_TO'"
    ).fetchone()[0]
    if old > 0:
        con.execute("DELETE FROM edges WHERE edge_type='SIMILAR_TO'")
        print(f"  Removed {old} old SIMILAR_TO edges")

    # ------------------------------------------------------------------ #
    # Insert new SIMILAR_TO edges
    # ------------------------------------------------------------------ #
    next_id = con.execute(
        "SELECT COALESCE(MAX(id), 0) FROM edges"
    ).fetchone()[0] + 1

    records = []
    for _, row in df.iterrows():
        records.append((
            next_id,
            "Structure", int(row["src_node"]),
            "Structure", int(row["tgt_node"]),
            "SIMILAR_TO",
            float(row["tmscore"]),   # weight = TM-score
            "Foldseek",
            0.8,                     # confidence
        ))
        next_id += 1

    con.executemany(
        "INSERT INTO edges VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        records,
    )
    print(f"  SIMILAR_TO edges inserted: {len(records):,}")

    # ------------------------------------------------------------------ #
    # Summary
    # ------------------------------------------------------------------ #
    total = con.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    by_type = con.execute(
        "SELECT edge_type, COUNT(*) FROM edges GROUP BY edge_type ORDER BY 2 DESC"
    ).fetchall()
    print()
    print("  Edge type breakdown:")
    for et, cnt in by_type:
        print(f"    {et}: {cnt:,}")
    print(f"  Total edges: {total:,}")

    con.close()
    print()
    print("Done. Next step: python3 scripts/materialize_graph.py \\")
    print("    --duckdb /home/.../atlas.duckdb \\")
    print("    --output-gpickle /home/.../atlas.gpickle")


if __name__ == "__main__":
    import os
    _base = Path(os.environ.get("PENSTACK_DATA", "/data"))
    ap = argparse.ArgumentParser()
    ap.add_argument("--duckdb",       type=Path, default=_base / "graphs/atlas.duckdb")
    ap.add_argument("--m8",           type=Path, default=_base / "foldseek/results.m8")
    ap.add_argument("--min-tm",       type=float, default=0.5)
    ap.add_argument("--max-per-node", type=int,   default=10)
    args = ap.parse_args()
    main(args.duckdb, args.m8, args.min_tm, args.max_per_node)
