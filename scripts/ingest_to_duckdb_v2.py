"""Clean ingestion pipeline v2 -- targets_v2.parquet, correct edge semantics."""
import duckdb, yaml, argparse, time, re
from pathlib import Path
import pandas as pd

DDL = """
CREATE TABLE IF NOT EXISTS nodes_system (
    id                   INTEGER PRIMARY KEY,
    name                 VARCHAR NOT NULL UNIQUE,
    type                 VARCHAR,
    subtype              VARCHAR,
    mechanism_bucket     VARCHAR CHECK (mechanism_bucket IN
        ('DSB_NUCLEASE','DSB_FREE_TRANSEST_RECOMBINASE','TRANSPOSASE','UNKNOWN')),
    rna_components       VARCHAR[],
    reference_doi        VARCHAR,
    notes                VARCHAR
);

CREATE TABLE IF NOT EXISTS nodes_protein (
    id                   INTEGER PRIMARY KEY,
    accession            VARCHAR NOT NULL UNIQUE,
    sequence             VARCHAR NOT NULL,
    length               INTEGER,
    organism_id          INTEGER,
    reviewed             BOOLEAN DEFAULT false,
    protein_name         VARCHAR
);

CREATE TABLE IF NOT EXISTS nodes_domain (
    id                   INTEGER PRIMARY KEY,
    accession            VARCHAR NOT NULL UNIQUE,
    name                 VARCHAR,
    source               VARCHAR DEFAULT 'Pfam',
    mechanism_bucket     VARCHAR
);

CREATE TABLE IF NOT EXISTS nodes_structure (
    id                   INTEGER PRIMARY KEY,
    accession            VARCHAR NOT NULL UNIQUE,
    source               VARCHAR,
    method               VARCHAR,
    resolution_A         DOUBLE,
    mean_plddt           DOUBLE
);

CREATE TABLE IF NOT EXISTS nodes_mechanism (
    id                   INTEGER PRIMARY KEY,
    name                 VARCHAR NOT NULL UNIQUE,
    bucket               VARCHAR,
    chemistry            VARCHAR,
    requires_host_repair BOOLEAN DEFAULT false
);

CREATE TABLE IF NOT EXISTS nodes_organism (
    id                   INTEGER PRIMARY KEY,
    ncbi_taxon_id        INTEGER NOT NULL UNIQUE,
    scientific_name      VARCHAR,
    lineage              VARCHAR[]
);

CREATE TABLE IF NOT EXISTS edges (
    id                   BIGINT PRIMARY KEY,
    source_type          VARCHAR NOT NULL,
    source_id            INTEGER NOT NULL,
    target_type          VARCHAR NOT NULL,
    target_id            INTEGER NOT NULL,
    edge_type            VARCHAR NOT NULL,
    weight               DOUBLE DEFAULT 1.0,
    provenance           VARCHAR,
    confidence           DOUBLE DEFAULT 1.0
);

CREATE INDEX IF NOT EXISTS idx_edges_source ON edges(source_type, source_id);
CREATE INDEX IF NOT EXISTS idx_edges_target ON edges(target_type, target_id);
CREATE INDEX IF NOT EXISTS idx_edges_type   ON edges(edge_type);
"""


def create_schema(con):
    con.execute(DDL)


def ingest_systems(con, yaml_path: Path):
    print("Ingesting systems...")
    data = yaml.safe_load(yaml_path.read_text())
    systems = data.get("systems", [])
    rows = []
    for i, s in enumerate(systems, 1):
        rna = s.get("rna_components", [])
        rows.append((
            i, s["name"], s.get("type"), s.get("subtype"),
            s.get("mechanism_bucket", "UNKNOWN"),
            rna if rna else None,
            s.get("reference_doi"), s.get("notes")
        ))
    con.executemany(
        "INSERT INTO nodes_system (id, name, type, subtype, mechanism_bucket, rna_components, reference_doi, notes) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
        rows
    )
    print(f"  Systems: {len(rows):,}")
    return {s["name"]: i for i, s in enumerate(systems, 1)}


def ingest_mechanisms(con):
    print("Ingesting mechanisms...")
    mechs = [
        (1, "DSB_NUCLEASE", "DSB_NUCLEASE", "metal-ion dependent hydrolysis", False),
        (2, "DSB_FREE_TRANSEST_RECOMBINASE", "DSB_FREE_TRANSEST_RECOMBINASE", "serine/tyrosine nucleophile", False),
        (3, "TRANSPOSASE", "TRANSPOSASE", "DDE catalysis", False),
    ]
    con.executemany(
        "INSERT INTO nodes_mechanism (id, name, bucket, chemistry, requires_host_repair) VALUES (?, ?, ?, ?, ?)",
        mechs
    )
    print(f"  Mechanisms: {len(mechs):,}")
    return {m[1]: m[0] for m in mechs}


def ingest_organisms(con, targets_parquet: Path):
    print("Ingesting organisms...")
    df = pd.read_parquet(targets_parquet, columns=["organism_id", "organism_name", "lineage_ids"])
    df = df.drop_duplicates(subset=["organism_id"]).reset_index(drop=True)
    
    def parse_lineage(lineage_str):
        if pd.isna(lineage_str) or str(lineage_str).strip() in ("", "None"):
            return None
        ids = re.findall(r"(\d+)", str(lineage_str))
        return ids if ids else None
    
    df["lineage"] = df["lineage_ids"].apply(parse_lineage)
    
    rows = []
    for i, row in df.iterrows():
        rows.append((
            i + 1, int(row["organism_id"]), str(row["organism_name"]),
            row["lineage"] if row["lineage"] else None
        ))
    
    con.executemany(
        "INSERT INTO nodes_organism (id, ncbi_taxon_id, scientific_name, lineage) VALUES (?, ?, ?, ?)",
        rows
    )
    print(f"  Organisms: {len(rows):,}")
    return {int(row["organism_id"]): i + 1 for i, row in df.iterrows()}


def ingest_proteins(con, targets_parquet: Path):
    print("Ingesting proteins...")
    df = pd.read_parquet(targets_parquet, columns=[
        "accession", "sequence", "length", "organism_id",
        "reviewed", "protein_name"
    ])
    
    rows = []
    for i, row in df.iterrows():
        rows.append((
            i + 1, str(row["accession"]), str(row["sequence"]),
            int(row["length"]) if pd.notna(row["length"]) else None,
            int(row["organism_id"]) if pd.notna(row["organism_id"]) else None,
            str(row["reviewed"]).lower() == "reviewed",
            str(row["protein_name"]) if pd.notna(row["protein_name"]) else None
        ))
    
    con.executemany(
        "INSERT INTO nodes_protein (id, accession, sequence, length, organism_id, reviewed, protein_name) VALUES (?, ?, ?, ?, ?, ?, ?)",
        rows
    )
    print(f"  Proteins: {len(rows):,}")
    return {str(row["accession"]): i + 1 for i, row in df.iterrows()}


def ingest_domains(con, whitelist_yaml: Path):
    print("Ingesting domains...")
    data = yaml.safe_load(whitelist_yaml.read_text())
    domains = data.get("domains", [])
    rows = []
    for i, d in enumerate(domains, 1):
        rows.append((
            i, str(d["accession"]), str(d.get("name", "")),
            "Pfam", str(d.get("mechanism_bucket", "UNKNOWN"))
        ))
    con.executemany(
        "INSERT INTO nodes_domain (id, accession, name, source, mechanism_bucket) VALUES (?, ?, ?, ?, ?)",
        rows
    )
    print(f"  Domains: {len(rows):,}")
    return {str(d["accession"]): i for i, d in enumerate(domains, 1)}


def ingest_structures(con, pdb_parquet: Path, afdb_parquet: Path):
    print("Ingesting structures...")
    
    pdb_df = pd.read_parquet(pdb_parquet)
    pdb_rows = []
    for i, row in pdb_df.iterrows():
        pdb_rows.append((
            i + 1, str(row.get("pdb_id", "")).upper(), "PDB",
            str(row.get("method", "")),
            float(row.get("resolution_A")) if pd.notna(row.get("resolution_A")) else None,
            None
        ))
    
    if pdb_rows:
        con.executemany(
            "INSERT INTO nodes_structure (id, accession, source, method, resolution_A, mean_plddt) VALUES (?, ?, ?, ?, ?, ?)",
            pdb_rows
        )
    
    max_id = len(pdb_rows)
    
    af_df = pd.read_parquet(afdb_parquet)
    af_rows = []
    for i, row in af_df.iterrows():
        af_rows.append((
            max_id + i + 1, f"AF-{str(row.get('accession', ''))}-F1", "AlphaFold",
            "Predicted", None,
            float(row.get("mean_plddt")) if pd.notna(row.get("mean_plddt")) else None
        ))
    
    if af_rows:
        con.executemany(
            "INSERT INTO nodes_structure (id, accession, source, method, resolution_A, mean_plddt) VALUES (?, ?, ?, ?, ?, ?)",
            af_rows
        )
    
    total = len(pdb_rows) + len(af_rows)
    print(f"  Structures: {total:,} (PDB: {len(pdb_rows):,}, AlphaFold: {len(af_rows):,})")
    return max_id + len(af_rows)


def build_edges_protein_domain(con, targets_parquet: Path, prot_map: dict, dom_map: dict, next_id: int):
    print("Building Protein->Domain edges...")
    df = pd.read_parquet(targets_parquet, columns=["accession", "xref_pfam"])
    df["xref_pfam"] = df["xref_pfam"].fillna("").astype(str)
    df = df[df["xref_pfam"] != ""]
    df["pfam_list"] = df["xref_pfam"].str.split(";")
    exploded = df.explode("pfam_list")
    exploded = exploded[exploded["pfam_list"].str.strip() != ""]
    exploded["pfam_acc"] = exploded["pfam_list"].str.strip()
    raw = exploded[["accession", "pfam_acc"]].drop_duplicates()
    
    raw["protein_id"] = raw["accession"].map(prot_map)
    raw["domain_id"] = raw["pfam_acc"].map(dom_map)
    raw = raw.dropna(subset=["protein_id", "domain_id"])
    
    n = len(raw)
    if n == 0:
        return 0, next_id
    
    edge_df = pd.DataFrame({
        "id": range(next_id, next_id + n),
        "source_type": "Protein",
        "source_id": raw["protein_id"].astype(int),
        "target_type": "Domain",
        "target_id": raw["domain_id"].astype(int),
        "edge_type": "HAS_DOMAIN",
        "weight": 1.0,
        "provenance": "UniProt/Pfam",
        "confidence": 1.0,
    })
    edge_df.to_parquet("/data/processed/edges_pd.parquet", index=False, compression="zstd")
    con.execute("INSERT INTO edges SELECT * FROM read_parquet('/data/processed/edges_pd.parquet')")
    print(f"  -> {n:,} edges")
    return n, next_id + n


def build_edges_system_protein(con, yaml_path: Path, sys_map: dict, prot_map: dict, next_id: int):
    print("Building System->Protein edges...")
    data = yaml.safe_load(yaml_path.read_text())
    systems = data.get("systems", [])
    
    rows = []
    for sys in systems:
        sid = sys_map.get(sys["name"])
        if not sid:
            continue
        for prot_acc in sys.get("proteins", []):
            pid = prot_map.get(prot_acc)
            if pid:
                rows.append((next_id, "System", sid, "Protein", pid, "HAS_PROTEIN", 1.0, "curated", 1.0))
                next_id += 1
    
    if rows:
        con.executemany(
            "INSERT INTO edges (id, source_type, source_id, target_type, target_id, edge_type, weight, provenance, confidence) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            rows
        )
    print(f"  -> {len(rows):,} edges")
    return len(rows), next_id


def build_edges_system_mechanism(con, yaml_path: Path, sys_map: dict, mech_map: dict, next_id: int):
    print("Building System->Mechanism edges...")
    data = yaml.safe_load(yaml_path.read_text())
    systems = data.get("systems", [])
    
    rows = []
    for sys in systems:
        sid = sys_map.get(sys["name"])
        bucket = sys.get("mechanism_bucket", "UNKNOWN")
        mid = mech_map.get(bucket)
        if sid and mid:
            rows.append((next_id, "System", sid, "Mechanism", mid, "USES_MECHANISM", 1.0, "curated", 1.0))
            next_id += 1
    
    if rows:
        con.executemany(
            "INSERT INTO edges (id, source_type, source_id, target_type, target_id, edge_type, weight, provenance, confidence) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            rows
        )
    print(f"  -> {len(rows):,} edges")
    return len(rows), next_id


def build_edges_structure_protein(con, targets_parquet: Path, prot_map: dict, next_id: int):
    print("Building Structure->Protein edges...")
    df = pd.read_parquet(targets_parquet, columns=["accession", "xref_pdb"])
    df = df[df["xref_pdb"].notna() & (df["xref_pdb"] != "") & (df["xref_pdb"] != "None")]
    df["xref_pdb"] = df["xref_pdb"].astype(str)
    
    structs = con.execute("SELECT id, accession FROM nodes_structure WHERE source = 'PDB'").fetchdf()
    struct_map = dict(zip(structs["accession"], structs["id"]))
    
    rows = []
    for _, row in df.iterrows():
        pid = prot_map.get(str(row["accession"]))
        if not pid:
            continue
        pdbs = [p.strip().upper() for p in str(row["xref_pdb"]).split(";") if p.strip()]
        for pdb in pdbs:
            sid = struct_map.get(pdb)
            if sid:
                rows.append((next_id, "Structure", int(sid), "Protein", pid, "STRUCTURE_OF", 1.0, "PDB", 1.0))
                next_id += 1
    
    if rows:
        con.executemany(
            "INSERT INTO edges (id, source_type, source_id, target_type, target_id, edge_type, weight, provenance, confidence) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
            rows
        )
    print(f"  -> {len(rows):,} edges")
    return len(rows), next_id


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data-dir", type=Path, required=True)
    p.add_argument("--pkg-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    
    print("=" * 60)
    print("GENOME-ATLAS Ingestion Pipeline v2")
    print("=" * 60)
    
    data_dir = args.data_dir
    pkg_dir = args.pkg_dir
    db_path = args.output
    
    db_path.parent.mkdir(parents=True, exist_ok=True)
    
    con = duckdb.connect(str(db_path))
    con.execute("SET memory_limit = '8GB'")
    
    print("")
    print("--- Creating schema ---")
    create_schema(con)
    
    print("")
    print("--- Ingesting nodes ---")
    sys_map = ingest_systems(con, pkg_dir / "genome_atlas" / "data" / "foundational_systems.yaml")
    mech_map = ingest_mechanisms(con)
    org_map = ingest_organisms(con, data_dir / "processed" / "targets_v2.parquet")
    prot_map = ingest_proteins(con, data_dir / "processed" / "targets_v2.parquet")
    dom_map = ingest_domains(con, pkg_dir / "genome_atlas" / "data" / "pfam_whitelist.yaml")
    n_structs = ingest_structures(con, data_dir / "processed" / "pdb_metadata.parquet",
                                   data_dir / "processed" / "alphafold_plddt.parquet")
    
    print("")
    print("--- Dropping indexes ---")
    con.execute("DROP INDEX IF EXISTS idx_edges_source")
    con.execute("DROP INDEX IF EXISTS idx_edges_target")
    con.execute("DROP INDEX IF EXISTS idx_edges_type")
    
    print("")
    print("--- Building edges ---")
    next_id = 1
    n_pd, next_id = build_edges_protein_domain(con, data_dir / "processed" / "targets_v2.parquet",
                                                prot_map, dom_map, next_id)
    n_sp, next_id = build_edges_system_protein(con, pkg_dir / "genome_atlas" / "data" / "foundational_systems.yaml",
                                                sys_map, prot_map, next_id)
    n_sm, next_id = build_edges_system_mechanism(con, pkg_dir / "genome_atlas" / "data" / "foundational_systems.yaml",
                                                  sys_map, mech_map, next_id)
    n_st, next_id = build_edges_structure_protein(con, data_dir / "processed" / "targets_v2.parquet",
                                                   prot_map, next_id)
    
    print("")
    print("--- Recreating indexes ---")
    con.execute("CREATE INDEX idx_edges_source ON edges(source_type, source_id)")
    con.execute("CREATE INDEX idx_edges_target ON edges(target_type, target_id)")
    con.execute("CREATE INDEX idx_edges_type ON edges(edge_type)")
    
    print("")
    print("=" * 60)
    print("FINAL COUNTS")
    print("=" * 60)
    for tbl in ["nodes_system", "nodes_protein", "nodes_domain", "nodes_structure", "nodes_mechanism", "nodes_organism"]:
        n = con.execute(f"SELECT COUNT(*) FROM {tbl}").fetchone()[0]
        print(f"  {tbl}: {n:,}")
    total_edges = con.execute("SELECT COUNT(*) FROM edges").fetchone()[0]
    print("")
    print(f"  Total edges: {total_edges:,}")
    for etype in ["HAS_DOMAIN", "HAS_PROTEIN", "USES_MECHANISM", "STRUCTURE_OF"]:
        n = con.execute(f"SELECT COUNT(*) FROM edges WHERE edge_type = '{etype}'").fetchone()[0]
        print(f"  {etype}: {n:,}")
    
    p = Path("/data/processed/edges_pd.parquet")
    if p.exists():
        p.unlink()
    
    con.execute("CHECKPOINT")
    con.close()
    print("")
    print("Done!")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
