import paramiko
import time

client = paramiko.SSHClient()
client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
client.connect('10.30.158.35', port=22, username='anees_22phd0670', password='78lh$K183#', timeout=30)

# Create a watcher script that runs extraction when download completes
watcher = """#!/bin/bash
LOG=~/pen-stack/data/afdb_pipeline.log
echo "$(date): AlphaFold pipeline watcher started" >> $LOG

# Wait for download tmux session to finish
while tmux has-session -t afdb-download 2>/dev/null; do
    sleep 30
done

echo "$(date): Download complete. Starting extraction..." >> $LOG

# Run extraction
bash ~/pen-stack/run_extraction.sh >> $LOG 2>&1

echo "$(date): Extraction complete." >> $LOG

# Cleanup: remove large tar files after extraction
# echo "$(date): Cleaning up tar files..." >> $LOG
# rm -f ~/pen-stack/data/raw/alphafold/*.tar
# echo "$(date): Cleanup complete." >> $LOG
"""

client.exec_command(f"cat > ~/pen-stack/watch_and_extract.sh <<'EOF'\n{watcher}\nEOF\nchmod +x ~/pen-stack/watch_and_extract.sh")

# Start watcher in tmux
client.exec_command('tmux kill-session -t afdb-watcher 2>/dev/null; sleep 1')
time.sleep(2)
client.exec_command("tmux new-session -d -s afdb-watcher 'bash ~/pen-stack/watch_and_extract.sh'")

print("Auto-pipeline set up:")
print("  - Download running in tmux: afdb-download")
print("  - Watcher running in tmux: afdb-watcher")
print("  - Watcher will auto-run extraction when download completes")

# Show current status
stdin, stdout, stderr = client.exec_command('tmux ls')
print(f"\nTmux sessions: {stdout.read().decode().strip()}")

client.close()
