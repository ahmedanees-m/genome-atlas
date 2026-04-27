import paramiko

client = paramiko.SSHClient()
client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
client.connect('10.30.158.35', port=22, username='anees_22phd0670', password='78lh$K183#', timeout=30)

print("=" * 70)
print("PEN-STACK VM PIPELINE STATUS")
print("=" * 70)

# Tmux sessions
stdin, stdout, stderr = client.exec_command('tmux ls 2>/dev/null')
sessions = stdout.read().decode().strip()
print(f"\nActive tmux sessions:\n{sessions or 'none'}")

# Step 6: UniProt
stdin, stdout, stderr = client.exec_command('ls ~/pen-stack/data/raw/uniprot/*.tsv.gz 2>/dev/null | wc -l')
uniprot_files = stdout.read().decode().strip()
stdin, stdout, stderr = client.exec_command('du -sh ~/pen-stack/data/raw/uniprot/')
uniprot_size = stdout.read().decode().strip()
stdin, stdout, stderr = client.exec_command('ls -lh ~/pen-stack/data/processed/targets_v1.parquet')
targets = stdout.read().decode().strip()
print(f"\n[Step 6] UniProt TSV Funnel")
print(f"  Files: {uniprot_files} | Size: {uniprot_size}")
print(f"  {targets}")

# Step 7: PDB
stdin, stdout, stderr = client.exec_command('ls ~/pen-stack/data/raw/pdb/*.pdb.gz 2>/dev/null | wc -l')
pdb_files = stdout.read().decode().strip()
stdin, stdout, stderr = client.exec_command('du -sh ~/pen-stack/data/raw/pdb/')
pdb_size = stdout.read().decode().strip()
stdin, stdout, stderr = client.exec_command('ls -lh ~/pen-stack/data/processed/pdb_metadata.parquet')
pdb_meta = stdout.read().decode().strip()
print(f"\n[Step 7] PDB Selective Download")
print(f"  Files: {pdb_files} | Size: {pdb_size}")
print(f"  {pdb_meta}")

# Step 8: AlphaFold
stdin, stdout, stderr = client.exec_command('ls ~/pen-stack/data/raw/alphafold/*.tar 2>/dev/null | wc -l')
af_files = stdout.read().decode().strip()
stdin, stdout, stderr = client.exec_command('du -sh ~/pen-stack/data/raw/alphafold/')
af_size = stdout.read().decode().strip()
stdin, stdout, stderr = client.exec_command('ls -lh ~/pen-stack/data/raw/alphafold/*.tar 2>/dev/null')
af_tars = stdout.read().decode().strip()
print(f"\n[Step 8] AlphaFold DB Download")
print(f"  Archives: {af_files} | Size: {af_size}")
for line in (af_tars or "").split('\n')[:5]:
    if line.strip():
        print(f"  {line}")
if af_tars and len(af_tars.split('\n')) > 5:
    print(f"  ... and {len(af_tars.split(chr(10))) - 5} more")

# Extraction status
stdin, stdout, stderr = client.exec_command('ls ~/pen-stack/data/raw/alphafold_targets/AF-*.pdb.gz 2>/dev/null | wc -l')
af_extracted = stdout.read().decode().strip()
stdin, stdout, stderr = client.exec_command('ls -lh ~/pen-stack/data/processed/alphafold_plddt.parquet 2>/dev/null')
af_plddt = stdout.read().decode().strip()
print(f"\n[Step 8] AlphaFold Extraction")
print(f"  Extracted structures: {af_extracted}")
print(f"  {af_plddt or 'pLDDT parquet: NOT YET CREATED'}")

# Pipeline log
stdin, stdout, stderr = client.exec_command('tail -5 ~/pen-stack/data/afdb_pipeline.log 2>/dev/null')
pipeline_log = stdout.read().decode().strip()
if pipeline_log:
    print(f"\nPipeline log tail:")
    for line in pipeline_log.split('\n'):
        print(f"  {line}")

print("\n" + "=" * 70)

client.close()
