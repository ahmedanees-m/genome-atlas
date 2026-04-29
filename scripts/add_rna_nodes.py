"""Add RNA nodes and System-HAS_RNA-RNA edges to the atlas DuckDB.

Creates nodes_rna table and populates it from the rnas: section of
foundational_systems.yaml, then wires HAS_RNA edges from each System
node to its constituent RNA nodes via the rna_components field.

Safe to re-run: uses ON CONFLICT DO NOTHING and skips existing edges.

Usage (on VM):
    python3 scripts/add_rna_nodes.py \
        --duckdb /home/anees_22phd0670/pen-stack/data/graphs/atlas.duckdb \
        --yaml   genome_atlas/data/foundational_systems.yaml
"""
import argparse
from pathlib import Path

import duckdb
import yaml


def main(duckdb_path: Path, yaml_path: Path):
    print("=" * 60)
    print("RNA Nodes Addition")
    print("=" * 60)

    cfg = yaml.safe_load(yaml_path.read_text())
    rnas = cfg.get("rnas", [])
    systems = cfg.get("systems", [])
    print(f"  YAML rnas:    {len(rnas)}")
    print(f"  YAML systems: {len(systems)}")

    con = duckdb.connect(str(duckdb_path))

    # ------------------------------------------------------------------ #
    # 1. Create nodes_rna table if absent
    # ------------------------------------------------------------------ #
    con.execute("""
        CREATE TABLE IF NOT EXISTS nodes_rna (
            id         INTEGER PRIMARY KEY,
            name       VARCHAR UNIQUE NOT NULL,
            rna_type   VARCHAR,
            length_nt  INTEGER,
            notes      VARCHAR
        )
    """)
    print("  nodes_rna table: ready")

    # ------------------------------------------------------------------ #
    # 2. Insert RNA nodes (idempotent via UNIQUE on name)
    # ------------------------------------------------------------------ #
    existing = {r[0] for r in con.execute("SELECT name FROM nodes_rna").fetchall()}
    next_id_row = con.execute("SELECT COALESCE(MAX(id), 0) FROM nodes_rna").fetchone()
    next_id = next_id_row[0] + 1

    inserted = 0
    for rna in rnas:
        name = rna["name"]
        if name in existing:
            print(f"    SKIP (exists): {name}")
            continue
        con.execute(
            "INSERT INTO nodes_rna VALUES (?, ?, ?, ?, ?)",
            [next_id, name, rna.get("rna_type"), rna.get("length_nt"), rna.get("notes")],
        )
        print(f"    INSERT RNA node: {name} (id={next_id})")
        next_id += 1
        inserted += 1
    print(f"  RNA nodes inserted: {inserted}")

    # ------------------------------------------------------------------ #
    # 3. Build rna_name → node_id lookup
    # ------------------------------------------------------------------ #
    rna_lookup = {
        r[0]: r[1]
        for r in con.execute("SELECT name, id FROM nodes_rna").fetchall()
    }

    # ------------------------------------------------------------------ #
    # 4. Build system_name → node_id lookup
    # ------------------------------------------------------------------ #
    sys_lookup = {
        r[0]: r[1]
        for r in con.execute("SELECT name, id FROM nodes_system").fetchall()
    }

    # ------------------------------------------------------------------ #
    # 5. Add HAS_RNA edges (System → RNA)
    # ------------------------------------------------------------------ #
    existing_edges = set(
        con.execute(
            "SELECT source_id, target_id FROM edges WHERE edge_type='HAS_RNA'"
        ).fetchall()
    )

    next_edge_id = con.execute("SELECT COALESCE(MAX(id), 0) FROM edges").fetchone()[0] + 1
    edges_added = 0
    edges_skipped = 0

    for sys in systems:
        sys_name = sys["name"]
        sys_id = sys_lookup.get(sys_name)
        if sys_id is None:
            print(f"    WARN: system '{sys_name}' not in DuckDB (not yet ingested)")
            continue

        for rna_name in sys.get("rna_components", []):
            rna_id = rna_lookup.get(rna_name)
            if rna_id is None:
                print(f"    WARN: RNA '{rna_name}' for system '{sys_name}' not in nodes_rna")
                continue

            if (sys_id, rna_id) in existing_edges:
                edges_skipped += 1
                continue

            con.execute(
                "INSERT INTO edges VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                [next_edge_id, "System", sys_id, "RNA", rna_id, "HAS_RNA",
                 1.0, "foundational_systems.yaml", 1.0],
            )
            existing_edges.add((sys_id, rna_id))
            next_edge_id += 1
            edges_added += 1

    print(f"  HAS_RNA edges added:   {edges_added}")
    print(f"  HAS_RNA edges skipped: {edges_skipped}")

    # ------------------------------------------------------------------ #
    # 6. Summary
    # ------------------------------------------------------------------ #
    rna_count  = con.execute("SELECT COUNT(*) FROM nodes_rna").fetchone()[0]
    edge_count = con.execute("SELECT COUNT(*) FROM edges WHERE edge_type='HAS_RNA'").fetchone()[0]
    total_edges = con.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    print()
    print("  Final state:")
    print(f"    nodes_rna rows : {rna_count}")
    print(f"    HAS_RNA edges  : {edge_count}")
    print(f"    Total edges    : {total_edges}")

    con.close()
    print()
    print("Done. Next step: python3 scripts/add_negative_controls.py")


if __name__ == "__main__":
    import os
    _base = Path(os.environ.get("PENSTACK_DATA", "/data"))
    ap = argparse.ArgumentParser()
    ap.add_argument("--duckdb", type=Path, default=_base / "graphs/atlas.duckdb")
    ap.add_argument("--yaml",   type=Path, default=Path("genome_atlas/data/foundational_systems.yaml"))
    args = ap.parse_args()
    main(args.duckdb, args.yaml)
