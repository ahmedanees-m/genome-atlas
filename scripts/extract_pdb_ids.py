"""Extract unique PDB IDs from targets table."""
import argparse
import sys
from pathlib import Path

import pandas as pd


def parse_pdb_ids(val):
    """Parse semicolon-separated PDB IDs from xref_pdb field."""
    if pd.isna(val) or not str(val).strip():
        return []
    ids = []
    for p in str(val).split(';'):
        p = p.strip().lower()
        if len(p) == 4 and p.isalnum():
            ids.append(p)
    return ids


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()

    df = pd.read_parquet(args.input)
    all_pdbs = set()
    for val in df['xref_pdb'].dropna():
        all_pdbs.update(parse_pdb_ids(val))

    sorted_ids = sorted(all_pdbs)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(sorted_ids))
    print(f"Wrote {len(sorted_ids):,} unique PDB IDs -> {args.output}")


if __name__ == "__main__":
    sys.exit(main())
