"""Extract metadata from downloaded PDB structure files."""
from __future__ import annotations

import argparse
import gzip
import re
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn


def parse_pdb_header(path: Path) -> dict | None:
    """Parse TITLE, RESOLUTION, EXPDTA, and HEADER date from a gzipped PDB file."""
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()
    except Exception:
        return None

    title_parts = []
    resolution = None
    method = None
    release_date = None

    for line in lines:
        record = line[:6].strip()

        if record == "TITLE":
            # TITLE     2 lines may be continued
            part = line[10:].strip()
            title_parts.append(part)

        elif record == "REMARK":
            if line[7:10] == "  2":
                # RESOLUTION
                m = re.search(r"RESOLUTION\.\s+([0-9.]+)\s+ANGSTROM", line)
                if m:
                    resolution = float(m.group(1))
            elif line[7:10] == "  1":
                # METHOD in REMARK 1
                if "RESOLUTION" not in line and method is None:
                    pass  # method usually in EXPDTA

        elif record == "EXPDTA":
            method = line[10:].strip()

        elif record == "HEADER":
            # HEADER    TRANSFERASE/TRANSPOSASE     18-OCT-23   8TZZ
            # The date is in columns 51-59
            date_str = line[50:59].strip()
            if date_str:
                try:
                    release_date = datetime.strptime(date_str, "%d-%b-%y").strftime("%Y-%m-%d")
                except ValueError:
                    pass

    title = " ".join(title_parts).strip()
    # Clean up continuation numbering like "2 STRUCTURE OF..."
    title = re.sub(r"^\d+\s+", "", title)

    return {
        "pdb_id": path.stem.replace(".pdb", "").lower(),
        "title": title,
        "resolution_A": resolution,
        "method": method,
        "release_date": release_date,
        "file_size_bytes": path.stat().st_size,
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--workers", type=int, default=2, help="Parallel workers (reduce if CPU-limited)")
    args = p.parse_args()

    files = sorted(args.input_dir.glob("*.pdb.gz"))
    if not files:
        print("ERROR: No .pdb.gz files found in", args.input_dir)
        return 1

    print(f"Extracting metadata from {len(files)} PDB files...")

    records = []
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
    ) as progress:
        task = progress.add_task("Parsing PDB headers", total=len(files))
        for fpath in files:
            meta = parse_pdb_header(fpath)
            if meta:
                records.append(meta)
            progress.advance(task)

    if not records:
        print("ERROR: No metadata extracted.")
        return 1

    df = pd.DataFrame(records)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.output, compression="zstd", index=False)

    print(f"\nWrote {len(df)} records to {args.output}")
    print(f"  With resolution: {df['resolution_A'].notna().sum()}")
    print(f"  With method: {df['method'].notna().sum()}")
    print(f"  With release_date: {df['release_date'].notna().sum()}")
    if df["method"].notna().any():
        print("  Methods:")
        for method, count in df["method"].value_counts().head(10).items():
            print(f"    {method}: {count}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
