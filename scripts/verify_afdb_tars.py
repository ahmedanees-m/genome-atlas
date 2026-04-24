"""Verify AlphaFold DB tar availability for top organisms."""
import pandas as pd
import requests
from pathlib import Path

df = pd.read_parquet("/data/processed/targets_v1.parquet")

# Top 30 organisms
org_counts = df.groupby(['organism_id', 'organism_name']).size().sort_values(ascending=False).head(30)

print("=== Verifying AlphaFold DB tar availability ===")
print("(Querying EBI FTP via HTTP)...")
print()

AFDB_BASE = "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/"

# Get directory listing
r = requests.get(AFDB_BASE, timeout=30)
if r.status_code != 200:
    print(f"ERROR: Cannot access AFDB FTP. HTTP {r.status_code}")
    exit(1)

# Parse HTML for tar files
import re
tar_files = re.findall(r'href="(UP[0-9]+_[0-9]+_[A-Za-z0-9_]+_v[0-9]+\.tar)"', r.text)
print(f"Total AFDB tars available: {len(tar_files)}")
print()

# Build lookup: taxid -> tar filename
tar_lookup = {}
for tar in tar_files:
    parts = tar.replace(".tar", "").split("_")
    if len(parts) >= 4:
        taxid = parts[2]  # UP000000625_83333_ECOLI_v6 -> 83333
        tar_lookup[taxid] = tar

print("Top organisms and AFDB coverage:")
print("-" * 100)
print(f"{'TaxID':>10} {'Count':>7} {'AF in targets':>13} {'Has Tar?':>8} {'Tar Filename':<50} {'Organism (truncated)'}")
print("-" * 100)

for (org_id, org_name), count in org_counts.items():
    org_id_str = str(org_id)
    af_count = df[(df.organism_id == org_id) & (df.xref_alphafolddb.notna())].shape[0]
    has_tar = org_id_str in tar_lookup
    tar_name = tar_lookup.get(org_id_str, "NOT FOUND")
    print(f"{org_id:>10} {count:>7,} {af_count:>13,} {'YES' if has_tar else 'NO':>8} {tar_name:<50} {org_name[:40]}")

print("-" * 100)

# Also check the execution plan's template organisms
print("\n=== Execution plan template organisms ===")
template = {
    "83333": "UP000000625_83333_ECOLI_v6",
    "1314": "UP000001215_1314_STRPY_v6",
    "1280": "UP000008816_93061_STAAU_v6",  # Note: key 1280 but filename has 93061
    "9606": "UP000005640_9606_HUMAN_v6",
}
for taxid, expected_tar in template.items():
    actual = tar_lookup.get(taxid, "NOT FOUND")
    count = df[df.organism_id == int(taxid)].shape[0] if taxid.isdigit() else 0
    status = "MATCH" if actual == expected_tar else "MISMATCH" if actual != "NOT FOUND" else "MISSING"
    print(f"  TaxID {taxid}: expected={expected_tar}, actual={actual}, status={status}, targets={count}")
