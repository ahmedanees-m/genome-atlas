"""Download AlphaFold DB archives for target organisms."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import requests
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn

AFDB_BASE = "https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/"

# Priority archives for genome-writing enzyme sources
# Covers: reviewed proteins + key bacterial model organisms
PRIORITY_ARCHIVES = [
    "swissprot_pdb_v6.tar",           # ALL reviewed proteins (28.6 GB)
    "UP000000625_83333_ECOLI_v6.tar",  # E. coli K-12 (0.5 GB)
    "UP000008816_93061_STAA8_v6.tar",  # S. aureus NCTC 8325 (0.3 GB)
    "UP000001014_99287_SALTY_v6.tar",  # Salmonella Typhimurium (0.5 GB)
    "UP000002438_208964_PSEAE_v6.tar", # P. aeruginosa (0.6 GB)
    "UP000000586_171101_STRR6_v6.tar", # S. pneumoniae (0.2 GB)
    "UP000007841_1125630_KLEPH_v6.tar", # K. pneumoniae (0.6 GB)
    "UP000001584_83332_MYCTU_v6.tar",  # M. tuberculosis (0.4 GB)
    "UP000000429_85962_HELPY_v6.tar",  # H. pylori (0.2 GB)
    "UP000000535_242231_NEIG1_v6.tar", # N. gonorrhoeae (0.2 GB)
]


def download_one(filename: str, out_dir: Path) -> bool:
    out = out_dir / filename
    if out.exists() and out.stat().st_size > 1000:
        print(f"  [skip] {filename} already exists ({out.stat().st_size / 1e9:.1f} GB)")
        return True

    url = AFDB_BASE + filename
    print(f"  [downloading] {filename} ...")
    try:
        r = requests.get(url, stream=True, timeout=60)
        r.raise_for_status()
        with open(out, "wb") as f:
            for chunk in r.iter_content(chunk_size=1048576):
                f.write(chunk)
        print(f"  [ok] {filename} -> {out.stat().st_size / 1e9:.1f} GB")
        return True
    except Exception as e:
        print(f"  [fail] {filename}: {e}")
        if out.exists():
            out.unlink()
        return False


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--archives", nargs="+", default=None, help="Specific archives to download")
    args = p.parse_args()

    args.output.mkdir(parents=True, exist_ok=True)
    archives = args.archives or PRIORITY_ARCHIVES

    print(f"Downloading {len(archives)} AlphaFold archives -> {args.output}")
    print(f"Expected total: ~{30:.1f} GB")
    print()

    success = 0
    fail = 0
    for filename in archives:
        if download_one(filename, args.output):
            success += 1
        else:
            fail += 1

    print(f"\nDone. Success: {success}, Failed: {fail}")
    total_gb = sum(f.stat().st_size for f in args.output.glob("*.tar")) / 1e9
    print(f"Total downloaded: {total_gb:.1f} GB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
