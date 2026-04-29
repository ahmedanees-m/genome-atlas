"""Analyze top organisms in targets_v1.parquet for AlphaFold coverage."""
import pandas as pd
from pathlib import Path
import requests

df = pd.read_parquet("/data/processed/targets_v1.parquet")

print("=== Top 30 Organisms by Protein Count ===")
org_counts = df.groupby(['organism_id', 'organism_name']).size().sort_values(ascending=False).head(30)
for (org_id, org_name), count in org_counts.items():
    reviewed = df[(df.organism_id == org_id) & (df.reviewed == 'reviewed')].shape[0]
    unreviewed = df[(df.organism_id == org_id) & (df.reviewed == 'unreviewed')].shape[0]
    with_pdb = df[(df.organism_id == org_id) & (df.xref_pdb.notna())].shape[0]
    with_af = df[(df.organism_id == org_id) & (df.xref_alphafolddb.notna())].shape[0]
    print(f"  {org_id:>10} | {count:>6,} total | {reviewed:>4} rev | {unreviewed:>6} unrev | {with_pdb:>4} PDB | {with_af:>6} AF | {org_name[:50]}")

print("\n=== Organisms with >1000 proteins ===")
big_orgs = org_counts[org_counts > 1000]
print(f"Count: {len(big_orgs)}")

print("\n=== PDB vs AlphaFold Coverage ===")
print(f"Total proteins: {len(df):,}")
print(f"With PDB: {df.xref_pdb.notna().sum():,} ({df.xref_pdb.notna().sum()/len(df)*100:.1f}%)")
print(f"With AlphaFold: {df.xref_alphafolddb.notna().sum():,} ({df.xref_alphafolddb.notna().sum()/len(df)*100:.1f}%)")
print(f"With neither: {(df.xref_pdb.isna() & df.xref_alphafolddb.isna()).sum():,}")
print(f"With both: {(df.xref_pdb.notna() & df.xref_alphafolddb.notna()).sum():,}")

print("\n=== AFDB Version Check ===")
r = requests.head("https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/", timeout=15, allow_redirects=True)
print(f"AFDB latest URL status: {r.status_code}")
if r.status_code == 200:
    print(f"Final URL: {r.url}")
