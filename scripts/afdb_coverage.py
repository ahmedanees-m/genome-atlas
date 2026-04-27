"""Calculate AFDB archive coverage for our targets."""
import pandas as pd
import requests
import json
import re

df = pd.read_parquet("/data/processed/targets_v1.parquet")

# Load AFDB metadata
r = requests.get('https://ftp.ebi.ac.uk/pub/databases/alphafold/download_metadata.json', timeout=30)
archives = r.json()

# Build taxid -> archive mapping
archive_by_taxid = {}
for item in archives:
    name = item.get('archive_name', '')
    parts = name.replace('.tar', '').split('_')
    if len(parts) >= 4 and parts[2].isdigit():
        taxid = parts[2]
        archive_by_taxid[taxid] = item

# For swissprot, use special key
for item in archives:
    if item.get('type') == 'swissprot' and 'pdb' in item.get('archive_name', ''):
        archive_by_taxid['swissprot_pdb'] = item

print("=== Archive Coverage Analysis ===\n")

# 1. SwissProt coverage (all reviewed proteins)
reviewed_targets = df[df.reviewed == 'reviewed']
print(f"SwissProt archive: {len(reviewed_targets):,} reviewed targets")

# 2. Per-archive coverage for proteome/global_health
results = []
for taxid, item in archive_by_taxid.items():
    if taxid == 'swissprot_pdb':
        continue
    try:
        targets = df[df.organism_id == int(taxid)].shape[0]
        if targets > 0:
            results.append((taxid, targets, item))
    except ValueError:
        pass

results.sort(key=lambda x: x[1], reverse=True)

print(f"\nMatching proteome/global_health archives: {len(results)}")
print(f"{'TaxID':>10} {'Targets':>8} {'Archive':<45} {'Species':<30} {'Size GB':>8}")
print("-" * 110)

total_targets = 0
total_size = 0
for taxid, targets, item in results:
    name = item.get('archive_name', 'N/A')
    species = item.get('species', 'N/A')[:28]
    size_gb = item.get('size_bytes', 0) / 1e9
    total_targets += targets
    total_size += size_gb
    print(f"{taxid:>10} {targets:>8,} {name:<45} {species:<30} {size_gb:>8.1f}")

print("-" * 110)
print(f"{'TOTAL':>10} {total_targets:>8,} {'':<45} {'':<30} {total_size:>8.1f}")

# 3. Best strategy
print(f"\n=== Recommended Download Strategy ===")
print(f"Option A: Download ALL matching archives ({len(results)} files, {total_size:.1f} GB)")
print(f"  -> Covers {total_targets:,} targets + {len(reviewed_targets):,} reviewed = {total_targets + len(reviewed_targets):,} total")

print(f"\nOption B: Download only top N archives + swissprot")
top5 = results[:5]
t5_targets = sum(r[1] for r in top5)
t5_size = sum(r[2].get('size_bytes', 0) for r in top5) / 1e9
print(f"  -> Top 5 archives ({t5_size:.1f} GB) + swissprot (28.6 GB) = {t5_size + 28.6:.1f} GB total")
print(f"  -> Covers {t5_targets:,} + {len(reviewed_targets):,} reviewed = {t5_targets + len(reviewed_targets):,} total")

# Show remaining targets without coverage
all_covered_taxids = set(r[0] for r in results) | {'swissprot_pdb'}
uncovered = df[~df.organism_id.astype(str).isin(all_covered_taxids)]
print(f"\nTargets NOT covered by any archive: {len(uncovered):,} ({len(uncovered)/len(df)*100:.1f}%)")
