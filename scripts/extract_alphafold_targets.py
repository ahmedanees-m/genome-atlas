"""Extract target proteins from AlphaFold DB tar archives."""
from __future__ import annotations

import argparse
import sys
import tarfile
from pathlib import Path

import pandas as pd
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--targets", type=Path, required=True)
    p.add_argument("--tar-dir", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    args = p.parse_args()

    # Load target accessions
    df = pd.read_parquet(args.targets)
    target_accs = set(df["accession"].astype(str).str.strip().unique())
    print(f"Total target accessions: {len(target_accs):,}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    tar_files = sorted(args.tar_dir.glob("*.tar"))
    if not tar_files:
        print("ERROR: No .tar files found in", args.tar_dir)
        return 1

    print(f"Scanning {len(tar_files)} tar archives...")

    total_extracted = 0
    total_matched = 0

    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
    ) as progress:
        task = progress.add_task("Extracting AlphaFold structures", total=len(tar_files))

        for tar_path in tar_files:
            progress.update(task, description=f"Scanning {tar_path.name}")
            try:
                with tarfile.open(tar_path, "r") as tar:
                    for member in tar.getmembers():
                        if not member.isfile():
                            continue
                        # AlphaFold filename: AF-<accession>-F1-model_v6.pdb.gz
                        name = member.name
                        if not name.startswith("AF-") or not name.endswith(".pdb.gz"):
                            continue
                        parts = name.replace(".pdb.gz", "").split("-")
                        if len(parts) < 2:
                            continue
                        acc = parts[1]
                        total_matched += 1
                        if acc in target_accs:
                            tar.extract(member, path=args.output_dir)
                            total_extracted += 1
            except Exception as e:
                print(f"WARNING: error processing {tar_path.name}: {e}")
            progress.advance(task)

    print(f"\nDone. Total PDB files in archives: {total_matched:,}")
    print(f"Extracted for targets: {total_extracted:,}")

    # Count actual extracted files
    extracted_files = list(args.output_dir.rglob("AF-*.pdb.gz"))
    print(f"Files in output dir: {len(extracted_files):,}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
