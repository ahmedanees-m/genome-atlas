"""Download PDB structures from RCSB via HTTP (Windows-compatible).

Replaces the rsync-based download_pdb_selective.sh from the execution plan,
since rsync is not available on Windows.
"""
from __future__ import annotations

import argparse
import gzip
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import requests
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn

RCSB_DOWNLOAD = "https://files.rcsb.org/download"


def download_one(pdb_id: str, out_dir: Path, fmt: str, max_retries: int = 3) -> tuple[str, bool]:
    """Download a single PDB structure. Returns (pdb_id, success)."""
    ext = "cif" if fmt == "cif" else "pdb"
    out = out_dir / f"{pdb_id}.{ext}.gz"
    if out.exists() and out.stat().st_size > 100:
        return pdb_id, True

    url = f"{RCSB_DOWNLOAD}/{pdb_id}.{ext}.gz"
    for attempt in range(max_retries):
        try:
            r = requests.get(url, timeout=30)
            if r.status_code == 404:
                return pdb_id, False  # Not found / deprecated
            r.raise_for_status()
            with open(out, "wb") as f:
                f.write(r.content)
            return pdb_id, True
        except Exception:
            if attempt < max_retries - 1:
                time.sleep(2 ** attempt)
    return pdb_id, False


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--ids", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--format", choices=["pdb", "cif"], default="pdb")
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--limit", type=int, default=0, help="Limit to first N IDs (0 = all)")
    args = p.parse_args()

    pdb_ids = [line.strip().lower() for line in args.ids.read_text().splitlines() if line.strip()]
    if args.limit > 0:
        pdb_ids = pdb_ids[:args.limit]
    args.output.mkdir(parents=True, exist_ok=True)

    print(f"Downloading {len(pdb_ids)} PDB structures ({args.format} format) -> {args.output}")
    print(f"Workers: {args.workers}")

    success = 0
    not_found = 0
    failed = 0

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
    ) as progress:
        task = progress.add_task("Downloading PDBs", total=len(pdb_ids))

        with ThreadPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(download_one, pid, args.output, args.format): pid for pid in pdb_ids}
            for future in as_completed(futures):
                pid, ok = future.result()
                if ok:
                    if (args.output / f"{pid}.{args.format}.gz").exists():
                        success += 1
                    else:
                        not_found += 1
                else:
                    not_found += 1
                progress.advance(task)

    print(f"\nDone. Success: {success}, Not found/deprecated: {not_found}, Failed: {failed}")
    print(f"Total downloaded: {len(list(args.output.glob('*.gz'))):,} files")
    total_mb = sum(f.stat().st_size for f in args.output.glob("*.gz")) / 1e6
    print(f"Total size: {total_mb:.1f} MB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
