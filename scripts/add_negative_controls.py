"""Add 500 negative-control proteins to the atlas DuckDB.

Selects 500 random unreviewed proteins from targets_v1 that are:
  - NOT in targets_v2 (not a genome-editing enzyme)
  - Have no Pfam accession overlap with the whitelist (filtered by xref_pfam IS NULL
    or containing no whitelisted domain prefix — simplified: no PF00 from whitelist)

Each negative control gets a USES_MECHANISM edge to a new NEGATIVE_CONTROL mechanism
node (id=4). Organisms not already in nodes_organism are inserted.

Safe to re-run: checks for accessions already in nodes_protein.

Usage (on VM):
    python3 scripts/add_negative_controls.py \
        --duckdb    /home/.../atlas.duckdb \
        --targets-v1 /home/.../processed/targets_v1.parquet \
        --targets-v2 /home/.../processed/targets_v2.parquet \
        --whitelist  genome_atlas/data/pfam_whitelist.yaml \
        --n 500 --seed 42
"""
import argparse
import random
from pathlib import Path

import duckdb
import yaml


def load_whitelist_accessions(yaml_path: Path) -> set:
    """Extract PF accessions from pfam_whitelist.yaml.

    Handles both the current layout (top-level 'domains' / 'auxiliary' lists
    of dicts with 'accession' keys) and older flat-list layouts.
    """
    cfg = yaml.safe_load(yaml_path.read_text())
    accs = set()
    # Current layout: 'domains' + 'auxiliary' sections
    for section in ("domains", "auxiliary", "primary_domains", "auxiliary_domains"):
        entries = cfg.get(section, [])
        if not isinstance(entries, list):
            continue
        for entry in entries:
            if isinstance(entry, dict) and "accession" in entry:
                accs.add(entry["accession"])
            elif isinstance(entry, str) and entry.startswith("PF"):
                accs.add(entry)
    return accs


def main(duckdb_path, t1_path, t2_path, whitelist_path, n, seed):
    print("=" * 60)
    print("Negative Controls Addition")
    print("=" * 60)

    whitelist = load_whitelist_accessions(whitelist_path)
    print(f"  Whitelist domains: {len(whitelist)}")

    con = duckdb.connect(str(duckdb_path))

    # ------------------------------------------------------------------ #
    # 1. Find already-present accessions
    # ------------------------------------------------------------------ #
    existing_accs = {
        r[0] for r in con.execute("SELECT accession FROM nodes_protein").fetchall()
    }
    print(f"  Existing nodes_protein: {len(existing_accs)}")

    # ------------------------------------------------------------------ #
    # 2. Sample negatives from targets_v1 \ targets_v2
    # ------------------------------------------------------------------ #
    print(f"  Sampling {n} negatives (seed={seed}) from targets_v1 \\ targets_v2 ...")
    t1 = str(t1_path)
    t2 = str(t2_path)

    # Build whitelist SQL filter: xref_pfam must not contain any whitelisted PF id
    # xref_pfam is a VARCHAR like "PF00589;PF00216" or NULL
    wl_filter_parts = [f"xref_pfam NOT LIKE '%{acc}%'" for acc in whitelist]
    wl_filter = " AND ".join(wl_filter_parts) if wl_filter_parts else "TRUE"

    # Fetch a large pool of unreviewed proteins with no Pfam hits, then filter in Python.
    # Avoid NOT IN (subquery) — it silently rejects all rows if subquery has NULLs.
    # Instead, fetch candidate pool and filter accessions in Python.
    query = f"""
        SELECT accession, sequence, length, organism_id, protein_name, organism_name
        FROM read_parquet('{t1}')
        WHERE reviewed = 'unreviewed'
          AND sequence IS NOT NULL
          AND length BETWEEN 50 AND 2000
          AND xref_pfam IS NULL
        LIMIT {n * 30}
    """
    candidates = con.execute(query).fetchdf()
    print(f"  Candidate pool (pre-filter): {len(candidates)}")

    # Python-side exclusion: remove proteins already in the atlas
    t2_accs = {r[0] for r in con.execute(
        f"SELECT accession FROM read_parquet('{t2}')"
    ).fetchall()}
    candidates = candidates[~candidates["accession"].isin(t2_accs)]
    candidates = candidates[~candidates["accession"].isin(existing_accs)]
    print(f"  After exclusion filter: {len(candidates)}")

    # Remove any already in nodes_protein
    candidates = candidates[~candidates["accession"].isin(existing_accs)]
    # Deterministic shuffle and take n
    rng = random.Random(seed)
    candidates = candidates.sample(n=min(n, len(candidates)), random_state=seed)
    print(f"  Final negatives selected: {len(candidates)}")

    # ------------------------------------------------------------------ #
    # 3. Ensure NEGATIVE_CONTROL mechanism node exists (id=4)
    # ------------------------------------------------------------------ #
    existing_mechs = {
        r[0] for r in con.execute("SELECT id FROM nodes_mechanism").fetchall()
    }
    if 4 not in existing_mechs:
        con.execute(
            "INSERT INTO nodes_mechanism VALUES (4, 'Negative Control', "
            "'NEGATIVE_CONTROL', 'none', false)"
        )
        print("  Added NEGATIVE_CONTROL mechanism (id=4)")
    else:
        print("  NEGATIVE_CONTROL mechanism already exists")

    # ------------------------------------------------------------------ #
    # 4. Ensure organisms exist in nodes_organism
    # ------------------------------------------------------------------ #
    existing_org_ids = {
        r[0] for r in con.execute("SELECT ncbi_taxon_id FROM nodes_organism").fetchall()
    }
    next_org_id = con.execute(
        "SELECT COALESCE(MAX(id), 0) FROM nodes_organism"
    ).fetchone()[0] + 1

    org_id_map = {}  # ncbi_taxon_id -> nodes_organism.id
    for row in con.execute(
        "SELECT id, ncbi_taxon_id FROM nodes_organism"
    ).fetchall():
        org_id_map[row[1]] = row[0]

    orgs_added = 0
    for _, cand_row in candidates.iterrows():
        taxon_id = int(cand_row["organism_id"]) if cand_row["organism_id"] else 0
        if taxon_id and taxon_id not in existing_org_ids:
            con.execute(
                "INSERT INTO nodes_organism VALUES (?, ?, ?, ?)",
                [next_org_id, taxon_id, str(cand_row["organism_name"]), []],
            )
            org_id_map[taxon_id] = next_org_id
            existing_org_ids.add(taxon_id)
            next_org_id += 1
            orgs_added += 1
    print(f"  Organisms added: {orgs_added}")

    # ------------------------------------------------------------------ #
    # 5. Insert negative control proteins into nodes_protein
    # ------------------------------------------------------------------ #
    next_prot_id = con.execute(
        "SELECT COALESCE(MAX(id), 0) FROM nodes_protein"
    ).fetchone()[0] + 1

    prots_added = 0
    neg_prot_ids = []

    for _, cand_row in candidates.iterrows():
        acc = cand_row["accession"]
        if acc in existing_accs:
            continue
        taxon_id = int(cand_row["organism_id"]) if cand_row["organism_id"] else 0
        # Map ncbi_taxon_id -> nodes_organism.id (default to 1 if missing)
        org_id = org_id_map.get(taxon_id, 1)

        con.execute(
            "INSERT INTO nodes_protein VALUES (?, ?, ?, ?, ?, ?, ?)",
            [
                next_prot_id,
                acc,
                str(cand_row["sequence"])[:65535],
                int(cand_row["length"]),
                org_id,
                False,   # reviewed=False (TrEMBL)
                str(cand_row["protein_name"])[:500] if cand_row["protein_name"] else "",
            ],
        )
        existing_accs.add(acc)
        neg_prot_ids.append(next_prot_id)
        next_prot_id += 1
        prots_added += 1

    print(f"  Protein nodes inserted: {prots_added}")

    # ------------------------------------------------------------------ #
    # 6. Add USES_MECHANISM edges for negatives → NEGATIVE_CONTROL (id=4)
    # ------------------------------------------------------------------ #
    next_edge_id = con.execute(
        "SELECT COALESCE(MAX(id), 0) FROM edges"
    ).fetchone()[0] + 1

    edges_added = 0
    for prot_id in neg_prot_ids:
        con.execute(
            "INSERT INTO edges VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            [next_edge_id, "Protein", prot_id, "Mechanism", 4,
             "USES_MECHANISM", 1.0, "negative_control", 0.5],
        )
        next_edge_id += 1
        edges_added += 1
    print(f"  USES_MECHANISM edges added: {edges_added}")

    # ------------------------------------------------------------------ #
    # 7. Summary
    # ------------------------------------------------------------------ #
    prot_total = con.execute("SELECT COUNT(*) FROM nodes_protein").fetchone()[0]
    edge_total = con.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    neg_edges  = con.execute(
        "SELECT COUNT(*) FROM edges WHERE edge_type='USES_MECHANISM' AND target_id=4"
    ).fetchone()[0]
    print()
    print("  Final state:")
    print(f"    nodes_protein total : {prot_total}")
    print(f"    USES_MECHANISM→NEG  : {neg_edges}")
    print(f"    Total edges         : {edge_total}")

    # ------------------------------------------------------------------ #
    # 8. Save negative control accessions for ESM-2 re-embedding
    # ------------------------------------------------------------------ #
    neg_accs = candidates["accession"].tolist()
    out_path = Path(str(duckdb_path)).parent / "negative_control_accessions.txt"
    out_path.write_text("\n".join(neg_accs) + "\n")
    print(f"    Accessions saved    : {out_path}")

    con.close()
    print()
    print("Done. Next step: bash scripts/run_foldseek.sh")


if __name__ == "__main__":
    # Default paths work both on VM directly and inside Docker
    # (Docker mounts ~/pen-stack/data -> /data)
    import os
    _base = Path(os.environ.get("PENSTACK_DATA", "/data"))
    ap = argparse.ArgumentParser()
    ap.add_argument("--duckdb",     type=Path,
                    default=_base / "graphs/atlas.duckdb")
    ap.add_argument("--targets-v1", type=Path,
                    default=_base / "processed/targets_v1.parquet")
    ap.add_argument("--targets-v2", type=Path,
                    default=_base / "processed/targets_v2.parquet")
    ap.add_argument("--whitelist",  type=Path,
                    default=Path("genome_atlas/data/pfam_whitelist.yaml"))
    ap.add_argument("--n",    type=int, default=500)
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    main(args.duckdb, args.targets_v1, args.targets_v2,
         args.whitelist, args.n, args.seed)
