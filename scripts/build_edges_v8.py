"""Build all edges v8 — drop indexes before bulk insert, recreate after."""
import duckdb, yaml, sys, time
from pathlib import Path
import pandas as pd

DB = Path("/data/graphs/atlas.duckdb")
PKG = Path("/pkg")
TARGETS = Path("/data/processed/targets_v1.parquet")

def get_next_id(con):
    r = con.execute("SELECT COALESCE(MAX(id), 0) FROM edges").fetchone()[0]
    return int(r) + 1 if r else 1

def build_protein_domain_edges(con, next_id: int):
    print("  Building protein-domain edges...")
    targets = pd.read_parquet(TARGETS, columns=["accession", "xref_pfam"])
    targets["xref_pfam"] = targets["xref_pfam"].fillna("").astype(str)
    targets = targets[targets["xref_pfam"] != ""]
    targets["pfam_list"] = targets["xref_pfam"].str.split(";")
    exploded = targets.explode("pfam_list")
    exploded = exploded[exploded["pfam_list"].str.strip() != ""]
    exploded["pfam_acc"] = exploded["pfam_list"].str.strip()
    raw = exploded[["accession", "pfam_acc"]].drop_duplicates()
    
    proteins = con.execute("SELECT id, accession FROM nodes_protein").fetchdf()
    domains = con.execute("SELECT id, accession FROM nodes_domain").fetchdf()
    
    prot_map = dict(zip(proteins["accession"], proteins["id"]))
    dom_map = dict(zip(domains["accession"], domains["id"]))
    
    raw["source_id"] = raw["accession"].map(prot_map)
    raw["target_id"] = raw["pfam_acc"].map(dom_map)
    raw = raw.dropna(subset=["source_id", "target_id"])
    
    n = len(raw)
    df = pd.DataFrame({
        "id": range(next_id, next_id + n),
        "source_type": "Protein",
        "source_id": raw["source_id"].astype(int),
        "target_type": "Domain",
        "target_id": raw["target_id"].astype(int),
        "edge_type": "HAS_DOMAIN",
        "weight": 1.0,
        "provenance": "UniProt/Pfam",
        "confidence": 1.0,
    })
    parquet_path = Path("/data/processed/edges_pd.parquet")
    df.to_parquet(parquet_path, index=False, compression="zstd")
    print(f"    -> {n:,} edges")
    return n, next_id + n

def build_system_protein_edges(con, next_id: int):
    print("  Building system-protein edges...")
    yaml_path = PKG / "genome_atlas" / "data" / "foundational_systems.yaml"
    data = yaml.safe_load(yaml_path.read_text())
    systems = data.get("systems", [])
    
    sys_nodes = con.execute("SELECT id, name FROM nodes_system").fetchdf()
    prot_nodes = con.execute("SELECT id, accession FROM nodes_protein").fetchdf()
    
    sys_map = dict(zip(sys_nodes["name"], sys_nodes["id"]))
    prot_map = dict(zip(prot_nodes["accession"], prot_nodes["id"]))
    
    rows = []
    for sys in systems:
        sid = sys_map.get(sys["name"])
        if not sid:
            continue
        for prot_acc in sys.get("proteins", []):
            pid = prot_map.get(prot_acc)
            if pid:
                rows.append((sid, pid))
    
    n = len(rows)
    df = pd.DataFrame({
        "id": range(next_id, next_id + n),
        "source_type": "System",
        "source_id": [r[0] for r in rows],
        "target_type": "Protein",
        "target_id": [r[1] for r in rows],
        "edge_type": "HAS_PROTEIN",
        "weight": 1.0,
        "provenance": "curated",
        "confidence": 1.0,
    })
    df.to_parquet("/data/processed/edges_sp.parquet", index=False, compression="zstd")
    print(f"    -> {n:,} edges")
    return n, next_id + n

def build_structure_protein_edges(con, next_id: int):
    print("  Building structure-protein edges...")
    proteins = pd.read_parquet(TARGETS, columns=["accession", "xref_pdb"])
    proteins = proteins[proteins["xref_pdb"].notna() & (proteins["xref_pdb"] != "")]
    proteins["xref_pdb"] = proteins["xref_pdb"].astype(str)
    
    structs = con.execute("SELECT id, accession FROM nodes_structure WHERE source = 'PDB'").fetchdf()
    struct_map = dict(zip(structs["accession"], structs["id"]))
    
    prot_nodes = con.execute("SELECT id, accession FROM nodes_protein").fetchdf()
    prot_map = dict(zip(prot_nodes["accession"], prot_nodes["id"]))
    
    rows = []
    for _, row in proteins.iterrows():
        pid = prot_map.get(row["accession"])
        if not pid:
            continue
        pdbs = [p.strip().upper() for p in row["xref_pdb"].split(";") if p.strip()]
        for pdb in pdbs:
            sid = struct_map.get(pdb)
            if sid:
                rows.append((int(sid), pid))
    
    n = len(rows)
    df = pd.DataFrame({
        "id": range(next_id, next_id + n),
        "source_type": "Structure",
        "source_id": [r[0] for r in rows],
        "target_type": "Protein",
        "target_id": [r[1] for r in rows],
        "edge_type": "STRUCTURE_OF",
        "weight": 1.0,
        "provenance": "PDB",
        "confidence": 1.0,
    })
    df.to_parquet("/data/processed/edges_st.parquet", index=False, compression="zstd")
    print(f"    -> {n:,} edges")
    return n, next_id + n

def build_mechanism_edges(con, next_id: int):
    print("  Building mechanism edges...")
    prots = pd.read_parquet(TARGETS, columns=["accession", "primary_mechanism_bucket"])
    prots = prots[prots["primary_mechanism_bucket"].notna() & (prots["primary_mechanism_bucket"] != "")]
    
    mechs = con.execute("SELECT id, lower(replace(trim(name), ' ', '_')) AS key FROM nodes_mechanism").fetchdf()
    mech_map = dict(zip(mechs["key"], mechs["id"]))
    
    prot_nodes = con.execute("SELECT id, accession FROM nodes_protein").fetchdf()
    prot_map = dict(zip(prot_nodes["accession"], prot_nodes["id"]))
    
    rows = []
    for _, row in prots.iterrows():
        pid = prot_map.get(row["accession"])
        if not pid:
            continue
        bucket = str(row["primary_mechanism_bucket"]).lower().replace(" ", "_")
        mid = mech_map.get(bucket)
        if mid:
            rows.append((pid, int(mid)))
    
    n = len(rows)
    df = pd.DataFrame({
        "id": range(next_id, next_id + n),
        "source_type": "Protein",
        "source_id": [r[0] for r in rows],
        "target_type": "Mechanism",
        "target_id": [r[1] for r in rows],
        "edge_type": "USES_MECHANISM",
        "weight": 1.0,
        "provenance": "curated",
        "confidence": 1.0,
    })
    df.to_parquet("/data/processed/edges_mech.parquet", index=False, compression="zstd")
    print(f"    -> {n:,} edges")
    return n, next_id + n

def build_organism_protein_edges(con, next_id: int):
    print("  Building organism-protein edges...")
    prots = con.execute("SELECT id, organism_id FROM nodes_protein WHERE organism_id IS NOT NULL").fetchdf()
    orgs = con.execute("SELECT id, ncbi_taxon_id FROM nodes_organism").fetchdf()
    
    org_map = dict(zip(orgs["ncbi_taxon_id"].astype(int), orgs["id"].astype(int)))
    rows = []
    for _, row in prots.iterrows():
        pid = int(row["id"])
        oid = org_map.get(int(row["organism_id"]))
        if oid:
            rows.append((oid, pid))
    
    n = len(rows)
    df = pd.DataFrame({
        "id": range(next_id, next_id + n),
        "source_type": "Organism",
        "source_id": [r[0] for r in rows],
        "target_type": "Protein",
        "target_id": [r[1] for r in rows],
        "edge_type": "HAS_PROTEIN",
        "weight": 1.0,
        "provenance": "UniProt",
        "confidence": 1.0,
    })
    df.to_parquet("/data/processed/edges_org.parquet", index=False, compression="zstd")
    print(f"    -> {n:,} edges")
    return n, next_id + n

def main():
    print("="*60)
    print("GENOME-ATLAS Edge Builder v8")
    print("="*60)
    
    if not DB.exists():
        print(f"ERROR: {DB} not found")
        sys.exit(1)
    
    con = duckdb.connect(str(DB))
    con.execute("SET memory_limit = '8GB'")
    next_id = get_next_id(con)
    print(f"Starting edge IDs from {next_id}")
    
    # Build all parquets
    print("\n--- Building parquet files ---")
    t0 = time.time()
    n_pd, next_id = build_protein_domain_edges(con, next_id)
    n_sp, next_id = build_system_protein_edges(con, next_id)
    n_st, next_id = build_structure_protein_edges(con, next_id)
    n_mech, next_id = build_mechanism_edges(con, next_id)
    n_org, next_id = build_organism_protein_edges(con, next_id)
    print(f"\nParquet build time: {time.time()-t0:.1f}s")
    
    # Drop indexes for faster insert
    print("\n--- Dropping indexes ---")
    con.execute("DROP INDEX IF EXISTS idx_edges_source")
    con.execute("DROP INDEX IF EXISTS idx_edges_target")
    con.execute("DROP INDEX IF EXISTS idx_edges_type")
    print("  Indexes dropped")
    
    # Insert all
    print("\n--- Bulk inserting into DuckDB ---")
    t0 = time.time()
    con.execute("INSERT INTO edges SELECT * FROM read_parquet('/data/processed/edges_pd.parquet')")
    print(f"  protein-domain: {time.time()-t0:.1f}s")
    t0 = time.time()
    con.execute("INSERT INTO edges SELECT * FROM read_parquet('/data/processed/edges_sp.parquet')")
    print(f"  system-protein: {time.time()-t0:.1f}s")
    t0 = time.time()
    con.execute("INSERT INTO edges SELECT * FROM read_parquet('/data/processed/edges_st.parquet')")
    print(f"  structure-protein: {time.time()-t0:.1f}s")
    t0 = time.time()
    con.execute("INSERT INTO edges SELECT * FROM read_parquet('/data/processed/edges_mech.parquet')")
    print(f"  mechanism: {time.time()-t0:.1f}s")
    t0 = time.time()
    con.execute("INSERT INTO edges SELECT * FROM read_parquet('/data/processed/edges_org.parquet')")
    print(f"  organism-protein: {time.time()-t0:.1f}s")
    
    # Recreate indexes
    print("\n--- Recreating indexes ---")
    t0 = time.time()
    con.execute("CREATE INDEX idx_edges_source ON edges(source_type, source_id)")
    con.execute("CREATE INDEX idx_edges_target ON edges(target_type, target_id)")
    con.execute("CREATE INDEX idx_edges_type ON edges(edge_type)")
    print(f"  Indexes recreated in {time.time()-t0:.1f}s")
    
    # Stats
    print("\n" + "="*60)
    print("FINAL EDGE COUNTS")
    print("="*60)
    for etype in ["HAS_DOMAIN", "HAS_PROTEIN", "STRUCTURE_OF", "USES_MECHANISM"]:
        cnt = con.execute(f"SELECT COUNT(*) FROM edges WHERE edge_type = '{etype}'").fetchone()[0]
        print(f"  {etype}: {cnt:,}")
    total = con.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    print(f"\n  TOTAL EDGES: {total:,}")
    
    # Cleanup
    for f in ["edges_pd", "edges_sp", "edges_st", "edges_mech", "edges_org"]:
        p = Path(f"/data/processed/{f}.parquet")
        if p.exists():
            p.unlink()
    
    con.execute("CHECKPOINT")
    con.close()
    print("\nDone!")
    return 0

if __name__ == "__main__":
    sys.exit(main())

