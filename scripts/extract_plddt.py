"""Extract pLDDT statistics from AlphaFold PDB structures."""
from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path

import pandas as pd
from Bio.PDB import PDBParser
from rich.progress import Progress, SpinnerColumn, BarColumn, TextColumn, TimeElapsedColumn


def parse_alphafold_plddt(path: Path) -> dict | None:
    """Parse pLDDT from CA atom B-factors in an AlphaFold PDB file."""
    try:
        parser = PDBParser(QUIET=True)
        if path.suffix == ".gz":
            with gzip.open(path, "rt") as f:
                structure = parser.get_structure("af", f)
        else:
            structure = parser.get_structure("af", str(path))
    except Exception:
        return None

    plddts = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    plddts.append(residue["CA"].get_bfactor())

    if not plddts:
        return None

    n = len(plddts)
    mean_p = sum(plddts) / n
    min_p = min(plddts)
    high_conf = sum(1 for p in plddts if p >= 70) / n

    # Extract accession from filename: AF-<acc>-F1-model_v6.pdb.gz
    acc = path.stem.replace(".pdb", "").split("-")[1] if "-" in path.stem else path.stem

    return {
        "accession": acc,
        "mean_plddt": round(mean_p, 2),
        "min_plddt": round(min_p, 2),
        "frac_high_confidence": round(high_conf, 4),
        "n_residues": n,
        "alphafold_path": str(path.name),
    }


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()

    files = list(args.input_dir.rglob("AF-*.pdb.gz"))
    if not files:
        print("ERROR: No AF-*.pdb.gz files found in", args.input_dir)
        return 1

    print(f"Computing pLDDT for {len(files)} AlphaFold structures...")

    records = []
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TextColumn("{task.completed}/{task.total}"),
        TimeElapsedColumn(),
    ) as progress:
        task = progress.add_task("Parsing pLDDT", total=len(files))
        for fpath in files:
            meta = parse_alphafold_plddt(fpath)
            if meta:
                records.append(meta)
            progress.advance(task)

    if not records:
        print("ERROR: No pLDDT records extracted.")
        return 1

    df = pd.DataFrame(records)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(args.output, compression="zstd", index=False)

    print(f"\nWrote {len(df)} records to {args.output}")
    print(f"  Mean pLDDT: {df['mean_plddt'].mean():.1f}")
    print(f"  Min pLDDT: {df['min_plddt'].min():.1f}")
    print(f"  High-conf fraction: {df['frac_high_confidence'].mean():.2%}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
